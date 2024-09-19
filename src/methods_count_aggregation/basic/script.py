import spatialdata as sd
import txsim as tx

## VIASH START
par = {
  'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr',
  'output': 'raw_counts.h5ad',
}
meta = {
  'name': 'basic'
}
## VIASH END

sdata = sd.read_zarr(par['input'])
df = sdata['transcripts'].compute() # TODO: Could optimize tx.preprocessing.generate_adata to work on spatialdata

adata = tx.preprocessing.generate_adata(df, cell_id_col='cell_id', gene_col='feature_name') #TODO: x and y refers to a specific coordinate system. Decide which space we want to use here. (probably should be handled in the previous assignment step)

# currently the function also saves the transcripts in the adata object, but this is not necessary here
del adata.uns['spots']
del adata.uns['pct_noise']

adata.write_h5ad(par['output'], compression="gzip")
