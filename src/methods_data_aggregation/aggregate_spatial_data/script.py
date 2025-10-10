import anndata as ad
import spatialdata as sd

## VIASH START
par = {
  'input_raw_sp': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
  'input_transcript_assignments': 'resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr',
  'input_qc_col': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_qc_col.h5ad',
  'input_spatial_corrected_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_corrected_counts.h5ad',
  'output': 'spatial_processed_complete.zarr',
}
meta = {
  'name': 'aggregate_spatial_data'
}
## VIASH END

#############
# Load data #
#############
sdata = sd.read_zarr(par['input_raw_sp'])
sdata_transcripts = sd.read_zarr(par['input_transcript_assignments'])
adata = ad.read_h5ad(par['input_spatial_corrected_counts'])
adata_qc_col = ad.read_h5ad(par['input_qc_col'])

###############
# Reduce data #
###############
# sdata
for key in list(sdata.labels.keys()):
  del sdata.labels[key]

for key in list(sdata.shapes.keys()):
  del sdata.shapes[key]

for key in list(sdata.points.keys()):
  del sdata.points[key]

for key in list(sdata.tables.keys()):
  if key != 'metadata':
    del sdata.tables[key]

# sdata_transcripts
for col in list(sdata_transcripts["transcripts"].columns):
  if col not in ['x', 'y', 'z', 'feature_name', 'cell_id', 'transcript_id']:
    del sdata_transcripts["transcripts"][col]

# adata
for col in list(adata.obs.columns):
  if col not in ['cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'n_counts', 'n_genes', 'volume', 'cell_type']:
    del adata.obs[col]

for col in list(adata.var.columns):
  if col not in ['gene_name', 'n_counts', 'n_cells']:
    del adata.var[col]

for layer in list(adata.layers.keys()):
  if layer not in ['counts', 'normalized', 'normalized_uncorrected']:
    del adata.layers[layer]

##############
# Assertions #
##############
assert len(sdata_transcripts["transcripts"]["cell_id"].unique()) in [adata.n_obs, adata.n_obs + 1], "Number of cells in transcripts and adata do not match"

##################
# Aggregate data #
##################
sdata["transcripts"] = sdata_transcripts["transcripts"]
adata.obs['passed_QC'] = adata_qc_col.obs['passed_QC']
sdata['counts'] = adata

#################
# Write output #
#################
sdata.write(par['output'])

