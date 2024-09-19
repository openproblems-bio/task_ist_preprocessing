import numpy as np
import anndata as ad
import spatialdata as sd
import os
import shutil

### VIASH START
par = {
  "input_scrnaseq": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_ist": "resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output_scrnaseq": "resources_test/task_ist_preprocessing/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_ist": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset.zarr"
}
### VIASH END

# Load the single-cell data
adata = ad.read_h5ad(par["input_scrnaseq"])

# Load the spatial data
sdata = sd.read_zarr(par["input_ist"])

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

# # NOTE: Arbitrary columns can lead to issues for anndata in R (specifically in the method script for 
# #       scrattch.mapping). In the future we'll need more columns (region matching filter) --> TODO
# # For now, only kept as a comment, should be fixed at the scrattch.mapping level
# adata.obs = adata.obs[["celltype"]]

# Filter spatial genes: #TODO: Add this here or as a requirement for raw spatialdata
# genes_sp = [g for g in genes_sp if not g.startswith("BLANK")]
# genes_sp = [g for g in genes_sp if not g.startswith("NegControl")]
# ...filter transcripts tables

# Save the single-cell data
adata.write_h5ad(par["output_scrnaseq"])

# remove directory if it exists
if os.path.exists(par["output_ist"]):
    shutil.rmtree(par["output_ist"])

# Save the spatial data
sdata.write(par["output_ist"], overwrite=True)
