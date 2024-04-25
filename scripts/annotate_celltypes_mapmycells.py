
from pathlib import Path
import argparse
import yaml
import anndata as ad
import json
import subprocess
import os
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate cells with celltype labels using tangram')

    parser.add_argument('-s', '--spatial', required=True, type=str,
	    help='Path to the spatial data h5ad')    
    parser.add_argument('-d', '--dissociated', required=True, type=str,
	    help='Path to the scRNAseq data h5ad')
    parser.add_argument('-t', '--temp_dir', required=True, type=str,
	    help='Temporary output directory')
    parser.add_argument('-o', '--output', required=True, type=str, 
        help='Output path of annotation csv file')    
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of method-specific parameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters')
    
    
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
    
    # Load default hyperparameter
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r',encoding='utf-8') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults["mapmycells"] 
        gparams_defaults = defaults["annotation_params"]
    
    # Parse parameters
    hyperparams = eval(args.hyperparams) if args.hyperparams else {}
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    groupparams = eval(args.groupparams) if args.groupparams else {}
    groupparams.update({k:v for k,v in gparams_defaults.items() if k not in groupparams})
    groupparams = {k:(v if v != "None" else None) for k,v in groupparams.items()}
    print("Hyperparameters:", hyperparams)
    print("Group parameters:", groupparams)
    
    # Read data
    adata_sc = ad.read_h5ad(args.dissociated)
    adata = ad.read_h5ad(args.spatial)
    
    
    #create taxonomy
    # hierarchy_levels = [hyperparams['class_level'], hyperparams['subclass_level'], hyperparams['cluster_level']]
    # hierarchy_string = json.dumps(hierarchy_levels)
    hierarchy_string = json.dumps([hyperparams['cell_type_column']])
    precomputed_stats_path = os.path.join(args.temp_dir,"precomputed_stats.h5ad")
    command = [
        "python",
        "-m",
        "cell_type_mapper.cli.precompute_stats_scrattch",
        "--h5ad_path",
        args.dissociated,  
        "--hierarchy",
        hierarchy_string,
        "--output_path",
        precomputed_stats_path
    ]

    subprocess.run(command)


    #create marker genes from spatial data genes
    genes= adata.var.index.values.tolist()
    data = {
        "None": genes
    }
    genes_file_path = os.path.join(args.temp_dir,"genes.json")
    with open(genes_file_path, "w") as json_file:
        json.dump(data, json_file, indent=2)

    #map cell types
    command = [
        "python",
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        args.spatial,  
        "--type_assignment.normalization",
        "log2CPM", #data is already normalized, if not set to raw
        
        "--precomputed_stats.path",
        precomputed_stats_path,
        "--query_markers.serialized_lookup",
        genes_file_path,
        "--csv_result_path",
        os.path.join(args.temp_dir,"results.csv"),
        "--extended_result_path",
        os.path.join(args.temp_dir,"extended_results.json"),
        #mapping onto a tree with only one taxonomic level; no bootstrapping
        "--flatten",
        "True",
        "--type_assignment.bootstrap_iteration", 
        "1",
        "--type_assignment.bootstrap_factor",
        "1.0"
 
    ]
    
    subprocess.run(command)

    #reformat output file
    annotation_df = pd.read_csv(os.path.join(args.temp_dir,"results.csv"), skiprows=3)
    annotation_df.drop([hyperparams['cell_type_column']+'_label', hyperparams['cell_type_column']+'_alias'], axis=1, inplace=True)
    annotation_df.rename(columns={hyperparams['cell_type_column']+'_name': 'celltype', hyperparams['cell_type_column']+'_correlation_coefficient': 'score'}, inplace=True)

    
    
    # Save annotation
    annotation_df.to_csv(args.output, index=False)