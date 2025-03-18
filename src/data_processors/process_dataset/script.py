import numpy as np
import anndata as ad
import spatialdata as sd
import os
import shutil

### VIASH START
par = {
  "input_sc": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_sp": "resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output_scrnaseq": "resources_test/task_ist_preprocessing/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_ist": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset.zarr"
}
### VIASH END

# Load the single-cell data
adata = ad.read_h5ad(par["input_sc"])

# Load the spatial data
sdata = sd.read_zarr(par["input_sp"])

# Subset the single-cell data to spatial genes
genes_sp = []
for key in sdata.tables.keys():
    # todo: var column names need to be updated to match the rest of openproblems
    genes_sp = genes_sp + sdata.tables[key].var_names.tolist()
genes_sp = list(np.unique(genes_sp))
adata = adata[:,adata.var["feature_name"].isin(genes_sp)].copy()

# Use feature names for adata instead of feature ids. convert to str
adata.var.reset_index(inplace=True, drop=True)
adata.var_names = adata.var["feature_name"].values.astype(str).tolist()

# store metadata to adata and sdata uns
metadata_uns_cols = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for col in metadata_uns_cols:
    orig_col = f"orig_{col}"
    if orig_col in adata.uns:
        adata.uns[orig_col] = adata.uns[col]
    adata.uns[col] = par[col]
    if orig_col in sdata.table.uns:
        sdata.table.uns[orig_col] = sdata.table.uns[col]
    sdata.table.uns[col] = par[col]

# Save the single-cell data
adata.write_h5ad(par["output_sc"], compression="gzip")

# remove directory if it exists
if os.path.exists(par["output_sp"]):
    shutil.rmtree(par["output_sp"])

# Save the spatial data
sdata.write(par["output_sp"], overwrite=True)
