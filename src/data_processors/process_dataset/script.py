import numpy as np
import scanpy as sc
import spatialdata as sd
import txsim as tx

### VIASH START
par = {
  "input_sc": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_sp": "resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
  "output_sc": "resources_test/task_ist_preprocessing/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_sp": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium/dataset.zarr"
}
### VIASH END

# Load the single-cell data
adata = sc.read(par["input_sc"])

# Load the spatial data
sdata = sd.read_zarr(par["input_sp"])

# Process single-cell data
adata = tx.preprocessing.normalize_sc(adata, layer="counts")
genes_sp = []
for key in sdata.points.keys():
    genes_sp = list(np.unique(genes_sp + sdata[key]["feature_name"].drop_duplicates().compute().tolist()))
adata = adata[:,adata.var["feature_name"].isin(genes_sp)]
assert len(adata.var["feature_name"]) == len(np.unique(adata.var["feature_name"])), "Gene symbols are not unique after subsetting to spatial genes"
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
