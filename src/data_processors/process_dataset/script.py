import numpy as np
import anndata as ad
import spatialdata as sd

### VIASH START
par = {
  "input_sc": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_sp": "resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output_sc": "resources_test/task_ist_preprocessing/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_sp": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset.zarr"
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
adata = adata[:,adata.var["feature_name"].isin(genes_sp)]

# Use feature names for adata instead of feature ids
adata.var_names = adata.var["feature_name"]

# # NOTE: Arbitrary columns can lead to issues for anndata in R (specifically in the method script for 
# #       scrattch.mapping). In the future we'll need more columns (region matching filter) --> TODO
# # For now, only kept as a comment, should be fixed at the scrattch.mapping level
# adata.obs = adata.obs[["celltype"]]

# Filter spatial genes: #TODO: Add this here or as a requirement for raw spatialdata
# genes_sp = [g for g in genes_sp if not g.startswith("BLANK")]
# genes_sp = [g for g in genes_sp if not g.startswith("NegControl")]
# ...filter transcripts tables

# Save the single-cell data
adata.write_h5ad(par["output_sc"])

# Save the spatial data
sdata.write(par["output_sp"], overwrite=True)
