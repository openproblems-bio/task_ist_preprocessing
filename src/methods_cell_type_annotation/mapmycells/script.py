import spatialdata as sd
import anndata as ad
import os
import subprocess
import json
## VIASH START
par = {
    'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
    'input_scrnaseq_reference': '.task_ist_preprocessing/resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/',
    "temp": './temp/mapmycells/',
    "output": 'output.h5ad' 
}
## VIASH END

#sdata = sd.read_zarr(par["input"])
if not os.path.exists(par["temp"]):
    os.makedirs(par["temp"])

adata = ad.read_h5ad(par['input_scrnaseq_reference'])
adata.X = adata.layers['counts']
adata.write_h5ad(f'{par["temp"]}sc_adata_processed.h5ad')


hierarchy_string = '["cell_type"]'

command = [
        "python",
        "-m",
        "cell_type_mapper.cli.precompute_stats_scrattch",
        "--h5ad_path",
        f'{par["temp"]}sc_adata_processed.h5ad',  
        "--hierarchy",
        hierarchy_string,
        "--output_path",
        f'{par["temp"]}precomputed_stats.h5ad'
    ]

subprocess.run(command)

genes= adata.var.index.values.tolist()
data = {
        "None": genes
}
genes_file_path = os.path.join(par["temp"],"genes.json")
with open(genes_file_path, "w") as json_file:
        json.dump(data, json_file, indent=2)

command = [
        "python",
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        #args.spatial,  
        "--type_assignment.normalization",
        "raw", 
        
        "--precomputed_stats.path",
        f'{par["temp"]}precomputed_stats.h5ad',
        "--query_markers.serialized_lookup",
        genes_file_path,
        "--csv_result_path",
        #os.path.join(args.temp_dir,"results.csv"),
        "--extended_result_path",
        #os.path.join(args.temp_dir,"extended_results.json"),
        #mapping onto a tree with only one taxonomic level; no bootstrapping
        "--flatten",
        "True",
        "--type_assignment.bootstrap_iteration", 
        "1",
        "--type_assignment.bootstrap_factor",
        "1.0"
 
    ]
    
#subprocess.run(command)

outputsdata = adata
outputsdata.write_h5ad(par['output'])