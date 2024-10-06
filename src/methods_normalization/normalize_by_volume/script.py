import anndata as ad
import txsim as tx

## VIASH START
par = {
  'input_spatial_aggregated_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad',
  'input_cell_volumes': 'resources_test/task_ist_preprocessing/mouse_brain_combined/cell_volumes.h5ad',
  'output': 'spatial_norm_counts.h5ad',
}
meta = {
  'name': 'normalize_by_volume',
}
## VIASH END

assert par['input_cell_volumes'] is not None, 'Volume input is required for this normalization method.'

print('Reading input files', flush=True)
adata = ad.read(par['input_spatial_aggregated_counts'])
adata_volume = ad.read(par['input_cell_volumes'])
assert adata.obs["cell_id"].astype(int).sort_values().equals(adata_volume.obs["cell_id"].astype(int).sort_values()), "Cell IDs do not match"
adata.obs["volume"] = adata_volume.obs.loc[adata.obs["cell_id"].astype(str), "volume"]

print('Normalizing by volume', flush=True)
 # TODO: Add scaling parameter, also as input to the script
tx.preprocessing.normalize_by_area(adata, area="volume")

adata.layers['normalized'] = adata.layers['lognorm'] #TODO: Check if openproblems "normalized" is understood as log normalized
del adata.layers['norm']
del adata.layers['lognorm']
#del adata.layers['raw'] #TODO: Probably better to add this, otherwise some downstream methods' tests will succeed by chance

print('Writing output', flush=True)
adata.write(par['output'])