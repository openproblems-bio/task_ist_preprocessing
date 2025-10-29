#!/usr/bin/env python3

# --- Fix for CuDNN version mismatch ------------------------------------------
# JAX wheels (e.g. jax[cuda12]) come with their own CUDA/cuDNN libraries.
# However, some systems set LD_LIBRARY_PATH to point to an older system-wide
# CUDA/cuDNN (e.g. 9.1). That can override JAX's bundled libraries and cause
# errors like:
#   "Loaded runtime CuDNN library: 9.1.0 but source was compiled with: 9.8.0"
#
# To prevent that, we remove LD_LIBRARY_PATH before importing JAX so it uses
# its own compatible, built-in CUDA/cuDNN stack.
# -----------------------------------------------------------------------------
import os
print("LD_LIBRARY_PATH before unset:", os.environ.get("LD_LIBRARY_PATH"), flush=True)
os.environ.pop("LD_LIBRARY_PATH", None)

# gpu check
import jax 
import jaxlib
print("GPU check", flush=True)
print("jax:", jax.__version__, flush=True)
print("jaxlib:", jaxlib.__version__, flush=True)
print("backend:", jax.default_backend(), flush=True)
print("devices:", jax.devices(), flush=True)
print("LD_LIBRARY_PATH:", os.environ.get("LD_LIBRARY_PATH"), flush=True)

import numpy as np
import anndata as ad
import scanpy as sc

import moscot as mt
from moscot.problems.space import MappingProblem

## VIASH START
par = {
   'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
   'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
   'output': 'spatial_with_celltypes.h5ad',
   'celltype_key': 'cell_type',
   'alpha': 0.8,
   'epsilon': 0.01,
   'tau': 1.0,
   'rank': 5000,
   'mapping_mode': 'max',
}
meta = {
   'name': 'moscot',
}
 ## VIASH END


# Optional parameter check: For this specific annotation method the par['input_spatial_normalized_counts'] and par['input_scrnaseq_reference'] are required
assert par['input_spatial_normalized_counts'] is not None, 'Spatial input is required for this annotation method.'
assert par['input_scrnaseq_reference'] is not None, 'Single cell input is required for this annotation method.'

# Read input
adata_sc = ad.read_h5ad(par['input_scrnaseq_reference'])
adata_sp = ad.read_h5ad(par['input_spatial_normalized_counts'])

# Adjust rank in case of small data set
if adata_sp.n_obs < 10000:
    print('Adjusting rank to -1 since data set is small (n_obs < 10k)', flush=True)
    par['rank'] = -1

# Check for normalized layer and centroid information
assert "normalized" in adata_sc.layers.keys(), 'Layer "normalized" is required for single-cell anndata'
assert "normalized" in adata_sp.layers.keys(), 'Layer "normalized" is required for spatial anndata'
assert "centroid_x" in adata_sp.obs and "centroid_y" in adata_sp.obs, '"Observation level columns "centroid_x" and "centroid_y" are required for spatial anndata'

# Use normalized layer and create spatial obsm
adata_sc.X = adata_sc.layers["normalized"]
adata_sp.X = adata_sp.layers["normalized"]
adata_sp.obsm["spatial"] = adata_sp.obs[["centroid_x", "centroid_y"]].to_numpy()

# Define mapping problem
mp = MappingProblem(adata_sc=adata_sc, adata_sp=adata_sp)

mp = mp.prepare(
    sc_attr={"attr": "layers", "key": "normalized"},
    xy_callback="local-pca",
)

mp = mp.solve(
    alpha=par['alpha'],
    epsilon=par['epsilon'],
    tau_a=par['tau'],
    tau_b=par['tau'],
    rank=par['rank'],
)

# Map annotations
anno_map_max = mp.annotation_mapping(
    mapping_mode=par['mapping_mode'],
    annotation_label=par['celltype_key'],
    source="src",
    forward=False,
)
adata_sp.obs[par['celltype_key']] = anno_map_max[par['celltype_key']].values

# Write output
adata_sp.write_h5ad(par['output'])
