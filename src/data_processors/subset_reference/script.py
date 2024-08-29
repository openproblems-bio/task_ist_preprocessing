import scanpy as sc
import numpy as np
import random

### VIASH START
par = {
    "input": "temp/datasets/allen_brain_cell_atlas/2023_Yao_mouse_brain_scRNAseq_10Xv2/tmp_dataset.h5ad",
    "output": "resources_test/common/2023_abca_Yao_mouse_brain_scRNAseq_10Xv2/dataset.h5ad",
    "n_cells": 500,
    "min_n_cells_per_cell_type": 50,
    "cell_type_key": "cell_type",
    "keep_cell_type_categories": None,
    "keep_features": None,
    "keep_tissue_categories": None,
    "seed": 0
}
### VIASH END

if par["seed"]:
    print(f">> Setting seed to {par['seed']}", flush=True)
    random.seed(par["seed"])
    np.random.seed(par["seed"])

print(">> Load data", flush=True)
adata = sc.read_h5ad(par["input"])

# Filtering by tissue categories if specified
if par.get("keep_tissue_categories"):
    print(f">> Filtering by tissue categories {par['keep_tissue_categories']}", flush=True)
    tissue_filt = adata.obs["tissue"].isin(par["keep_tissue_categories"])
    adata = adata[tissue_filt].copy()

# Filtering by cell type categories if specified
if par.get("keep_cell_type_categories"):
    print(f">> Filtering by cell type categories {par['keep_cell_type_categories']}", flush=True)
    cell_type_filt = adata.obs[par["cell_type_key"]].isin(par["keep_cell_type_categories"])
    adata = adata[cell_type_filt].copy()

# Subsampling cells based on the minimum number of cells per cell type
cell_types = adata.obs[par["cell_type_key"]].unique()
selected_indices = []

n_per_cell_type = {}
for ct in cell_types:
    cell_type_indices = adata.obs_names[adata.obs[par["cell_type_key"]] == ct]
    n_cells_in_type = len(cell_type_indices)
    n_per_cell_type[ct] = min(n_cells_in_type, par["min_n_cells_per_cell_type"])
    selected_indices.extend(np.random.choice(cell_type_indices, n_per_cell_type[ct], replace=False))

# Cap the number of cells to sample per cell type by the total number of cells to sample
#TODO: instead of random choice below, adjust that sum(n_per_cell_type.values()) <= par["n_cells"]

# Ensure we don't exceed the total number of cells to sample
selected_indices = np.random.choice(selected_indices, par["n_cells"], replace=False)
adata = adata[selected_indices].copy()

# Subsampling features (genes) if specified
if par.get("keep_features"):
    print(f">> Keeping specified features {par['keep_features']}", flush=True)
    feature_filt = adata.var_names.isin(par["keep_features"])
    adata = adata[:, feature_filt].copy()

# Remove empty observations and features
print(">> Removing empty observations and features", flush=True)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=2)

# Update dataset_id or add relevant metadata if needed
print(">> Update metadata", flush=True)
adata.uns["dataset_id"] = adata.uns.get("dataset_id", "dataset") + "_subset"

# Write the subsetted data to output file
print(">> Writing data", flush=True)
adata.write_h5ad(par["output"])

print(">> Subsetting complete", flush=True)
