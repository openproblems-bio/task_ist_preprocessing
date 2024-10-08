import anndata as ad
import txsim as tx

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad',
  'min_counts': 10,
  'min_cell_percentage': 0.8,
  'output': 'spatial_qc_col.h5ad',
}
meta = {
  'name': 'basic_qc_filter',
}
## VIASH END
    
# Read input
print('Reading input files', flush=True)
adata = ad.read_h5ad(par['input'])

# Apply gene efficiency correction
print('Get QC filter', flush=True)
tx.preprocessing.filter_cells(
    adata, min_counts = par['min_counts'], min_cell_percentage = par['min_cell_percentage'], obs_key = "passed_QC"
)
# NOTE: this function makes use of the column adata.obs['n_counts'] which is calculated in the count aggregation step. 
# TODO: Add n_counts column to according viash config

# Write output
print('Writing output', flush=True)
adata_obs_col = ad.AnnData(obs=adata.obs[["passed_QC"]])
adata_obs_col.write(par['output'])
