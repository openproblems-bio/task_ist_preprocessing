import anndata as ad
import txsim as tx

## VIASH START
par = {
  'input_spatial_aggregated_counts': 'spatial_raw_counts.h5ad',
  'input_cell_volumes': 'cell_volumes.h5ad',
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

print('Writing output', flush=True)
adata.write(par['output'])