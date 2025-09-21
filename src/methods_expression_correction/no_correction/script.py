import anndata as ad

## VIASH START
par = {
  'input_spatial_with_cell_types': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_with_celltypes.h5ad',
  'output': 'spatial_corrected.h5ad',
}
## VIASH END
    
# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input_spatial_with_cell_types'])
adata_sp.layers["normalized_uncorrected"] = adata_sp.layers["normalized"]

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])
