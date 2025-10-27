import spatialdata as sd
import anndata as ad
import os
import subprocess
import json
import pandas as pd
## VIASH START
par = {
    'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
    'input_scrnaseq_reference': './resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad',
    "temp": './task_ist_preprocessing/temp/mapmycells/',
    'cell_type_column': 'celltype',
    "output": 'output.h5ad'
}
## VIASH END


if not os.path.exists(par["temp"]):
    os.makedirs(par["temp"])

sdata= ad.read_h5ad(par['input_spatial_normalized_counts'])
adata = ad.read_h5ad(par['input_scrnaseq_reference'])

if "counts" in adata.layers:
    adata.X = adata.layers["counts"]

mapping = dict(zip(adata.var["feature_name"], adata.var["feature_id"]))
sdata.var.index = sdata.var.index.map(mapping)
adata.var_names_make_unique()  
sdata.var.index = sdata.var.index.astype(str)
sdata.var_names_make_unique()
common_genes = list(set(sdata.var.index).intersection(adata.var.index))

adata = adata[:, common_genes]



sc_path = f'{par["temp"]}sc_adata_processed.h5ad'
if os.path.exists(sc_path):
    os.remove(sc_path)
adata.write_h5ad(sc_path)

# 🧹 Check before writing sp_processed.h5ad
sp_path = f'{par["temp"]}sp_processed.h5ad'
if os.path.exists(sp_path):
    os.remove(sp_path)
sdata[:, common_genes].write_h5ad(sp_path)



hierarchy_string = '["cell_type"]'

precomputed_path = f'{par["temp"]}precomputed_stats.h5ad'
if os.path.exists(precomputed_path):
    os.remove(precomputed_path)

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



data = {"None": common_genes}
genes_file_path = os.path.join(par["temp"],"genes.json")
if os.path.exists(genes_file_path):
    os.remove(genes_file_path)

with open(genes_file_path, "w") as json_file:
        json.dump(data, json_file, indent=2)


command = [
        "python",
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        f'{par["temp"]}sp_processed.h5ad',  
        "--type_assignment.normalization",
        "log2CPM", 
        "--precomputed_stats.path",
        f'{par["temp"]}precomputed_stats.h5ad',
        "--query_markers.serialized_lookup",
        genes_file_path,
        "--csv_result_path",
        os.path.join(par["temp"],"results.csv"),
        "--extended_result_path",
        os.path.join(par["temp"], "extended_results.json"),
        #mapping onto a tree with only one taxonomic level; no bootstrapping
        "--flatten",
        "True",
        "--type_assignment.bootstrap_iteration", 
        "1",
        "--type_assignment.bootstrap_factor",
        "1.0"
 
    ]

results_csv = os.path.join(par["temp"], "results.csv")
if os.path.exists(results_csv):
    os.remove(results_csv)

extended_json = os.path.join(par["temp"], "extended_results.json")
if os.path.exists(extended_json):
    os.remove(extended_json)

    
subprocess.run(command)

label_name = 'cell_type'
annotation_df = pd.read_csv(f'{par["temp"]}results.csv', skiprows=3)

sdata.obs[par['cell_type_column']] = annotation_df['cell_type_name']

sdata.write_h5ad(par['output'])