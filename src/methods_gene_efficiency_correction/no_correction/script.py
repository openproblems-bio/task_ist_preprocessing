import anndata as ad

## VIASH START
par = {
  'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_corrected.h5ad',
  'output': 'spatial_gene_efficiency_corrected.h5ad',
}
## VIASH END

# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input'])
# Preserve the pre-correction normalized layer only if an upstream stage hasn't set it.
if "normalized_uncorrected" not in adata_sp.layers:
    adata_sp.layers["normalized_uncorrected"] = adata_sp.layers["normalized"]

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])
