import anndata as ad
import os
import subprocess
import json
import pandas as pd
from pathlib import Path 
## VIASH START
par = {
    'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
    'input_scrnaseq_reference': 'resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad',
    'cell_type_column': 'cell_type',
    "output": 'output.h5ad'
}
meta = { "temp_dir": './tmp/'}

## VIASH END

TMP_DIR = Path(meta["temp_dir"] or "/tmp/")
TMP_DIR.mkdir(parents=True, exist_ok=True)

sdata = ad.read_h5ad(par['input_spatial_normalized_counts'])
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
sc_path = os.path.join(meta["temp_dir"],"sc_adata_processed.h5ad")
adata.write_h5ad(sc_path)
sp_path = os.path.join(meta["temp_dir"],"sp_processed.h5ad")
sdata[:, common_genes].write_h5ad(sp_path)



precomputed_path = os.path.join(meta["temp_dir"],"precomputed_stats.h5ad")

command = [
        "python",
        "-m",
        "cell_type_mapper.cli.precompute_stats_scrattch",
        "--h5ad_path",
        sc_path,  
        "--hierarchy",
        "['cell_type']",
        "--output_path",
       precomputed_path
    ]

subprocess.run(command)

data = {"None": common_genes}
genes_file_path = os.path.join(meta["temp_dir"],"genes.json")
with open(genes_file_path, "w") as json_file:
        json.dump(data, json_file, indent=2)
command = [
        "python",
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        sp_path,  
        "--type_assignment.normalization",
        "log2CPM", 
        "--precomputed_stats.path",
        precomputed_path,
        "--query_markers.serialized_lookup",
        genes_file_path,
        "--csv_result_path",
        os.path.join(meta["temp_dir"],"results.csv"),
        "--extended_result_path",
        os.path.join(meta["temp_dir"], "extended_results.json"),
        "--flatten",
        "True",
        "--type_assignment.bootstrap_iteration", 
        "1",
        "--type_assignment.bootstrap_factor",
        "1.0"

    ]

subprocess.run(command)
annotation_df = pd.read_csv(os.path.join(meta["temp_dir"],"results.csv"), skiprows=3)
sdata.obs[par['cell_type_column']] = list(annotation_df['cell_type_label'])



# Loop through all files in the folder
for file_path in [
     sc_path,
     sp_path,
     precomputed_path,
     genes_file_path,
     os.path.join(meta["temp_dir"],"results.csv"),
     os.path.join(meta["temp_dir"], "extended_results.json")
]:
        if os.path.isfile(file_path):
            os.remove(file_path)


sdata.write_h5ad(par['output'])
