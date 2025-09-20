import anndata as ad
import scanpy as sc

## VIASH START
par = {
  'input_spatial_aggregated_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad',
  'output': 'spatial_norm_counts.h5ad',

  'target_sum': None,
}
## VIASH END

print('Reading input files', flush=True)
adata = ad.read(par['input_spatial_aggregated_counts'])

print('Normalizing by total counts', flush=True)
adata.layers['normalized'] = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, layer="normalized", target_sum=par["target_sum"])
sc.pp.log1p(adata, layer="normalized")

print('Writing output', flush=True)
adata.write(par['output'])