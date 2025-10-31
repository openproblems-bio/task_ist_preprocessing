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

# Subset single-cell and spatial data to shared genes
sp_genes = sdata['transcripts']['feature_name'].unique().compute().tolist()
sc_genes = adata.var["feature_name"].unique().tolist()
shared_genes = list(set(sp_genes) & set(sc_genes))
sdata['transcripts'] = sdata['transcripts'].loc[sdata['transcripts']['feature_name'].isin(shared_genes)]
adata = adata[:,adata.var["feature_name"].isin(shared_genes)].copy()

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

# Correct the feature_key attribute in sdata if needed
# NOTE: it would have been better to do this in the loader scripts, but this way the datasets don't need to be re-downloaded
if "feature_key" in sdata['transcripts'].attrs["spatialdata_attrs"]:
    feature_key = sdata['transcripts'].attrs["spatialdata_attrs"]["feature_key"]
    if feature_key != "feature_name":
        sdata['transcripts'].attrs["spatialdata_attrs"]["feature_key"] = "feature_name"

# Save the single-cell data
adata.write_h5ad(par["output_sc"], compression="gzip")

# remove directory if it exists
if os.path.exists(par["output_sp"]):
    shutil.rmtree(par["output_sp"])

# Save the spatial data
sdata.write(par["output_sp"], overwrite=True)
