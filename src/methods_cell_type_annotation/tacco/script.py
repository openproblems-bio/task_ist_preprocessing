#!/usr/bin/env python3

import anndata as ad
import numpy as np
import tacco

## VIASH START
par = {
  'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
  'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
  'output': 'spatial_with_celltypes.h5ad',
  'celltype_key': 'cell_type',
}
meta = {
  'name': 'tacco',
}
## VIASH END

# Optional parameter check: For this specific annotation method the par['input_spatial_normalized_counts'] and par['input_scrnaseq_reference'] are required
assert par['input_spatial_normalized_counts'] is not None, 'Spatial input is required for this annotation method.'
assert par['input_scrnaseq_reference'] is not None, 'Single cell input is required for this annotation method.'

# Read input
adata_sp = ad.read_h5ad(par['input_spatial_normalized_counts'])
adata_sc = ad.read_h5ad(par['input_scrnaseq_reference'])

# Switch to raw counts
adata_sp.X = adata_sp.layers['counts']
adata_sc.X = adata_sc.layers['counts']

# Run tacco
cell_type_assignment = tacco.tl.annotate(
    adata=adata_sp,
    reference=adata_sc,
    annotation_key=par['celltype_key']
)

# Tacco stores the cell type proportions in a n_obs x n_celltypes matrix, so we have to extract the celltype with highest consensus
cell_types = cell_type_assignment.columns
highest_score_idx = np.argmax(cell_type_assignment, axis=1)
adata_sp.obs[par['celltype_key']] = cell_types[highest_score_idx]

# Write output
adata_sp.write_h5ad(par['output'])
