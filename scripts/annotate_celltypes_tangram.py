
from pathlib import Path
import argparse
import yaml
import anndata as ad

#from _ctannotation import run_tangram 
from txsim.preprocessing._ctannotation import run_tangram




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate cells with celltype labels using tangram')

    parser.add_argument('-s', '--spatial', required=True, type=str,
	    help='Path to the spatial data h5ad')    
    parser.add_argument('-d', '--dissociated', required=True, type=str,
	    help='Path to the scRNAseq data h5ad')
    parser.add_argument('-o', '--output', required=True, type=str, 
        help='Output path of annotation csv file')    
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of method-specific parameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters')
    
    
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Load default hyperparameter
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r',encoding='utf-8') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults["tangram"] 
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
    
    # Run cell type annotation 
    annotation_df = run_tangram(adata, adata_sc, sc_ct_labels=hyperparams["sc_ct_labels"],  mode = hyperparams["mode"],num_epochs = hyperparams["num_epochs"],device = hyperparams["device"])


    
    
    # Save annotation
    annotation_df.obs.to_csv(args.output, index=False)