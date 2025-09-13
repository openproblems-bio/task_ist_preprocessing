import numpy as np
import anndata as ad

## VIASH START
par = {
    'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
    'output': 'output.h5ad',
}
## VIASH END

print('Read input_scrnaseq_reference', flush=True)
adata = ad.read_h5ad(par['input_scrnaseq_reference'])

print("Add required layers, obs and var columns for spatial data", flush=True)
adata.layers["normalized_uncorrected"] = adata.layers["normalized"]
adata.obs["cell_id"] = adata.obs.index
adata.obs["centroid_x"] = 0
adata.obs["centroid_y"] = 0
adata.obs["centroid_z"] = 0
adata.obs["n_counts"] = np.array(adata.layers["counts"].sum(axis=1))[:,0]
adata.obs["n_genes"] = np.array((adata.layers["counts"] > 0).sum(axis=1))[:,0]
adata.obs["volume"] = 1
adata.var["gene_name"] = adata.var.index
adata.var["n_counts"] = np.array(adata.layers["counts"].sum(axis=0))[0,:]
adata.var["n_cells"] = np.array((adata.layers["counts"] > 0).sum(axis=0))[0,:]

print("Store outputs", flush=True)
adata.write_h5ad(par['output'], compression='gzip')
