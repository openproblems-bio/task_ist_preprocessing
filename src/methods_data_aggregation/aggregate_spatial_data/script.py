import anndata as ad
import geopandas as gpd
import sopa
import spatialdata as sd
from shapely.geometry import MultiPoint
from spatialdata.models import ShapesModel
from sopa.utils import copy_transformations

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
  if key not in ['metadata', 'table']:
    del sdata.tables[key]

# raw_ist.zarr stores the metadata table as 'table'; rename to match the output spec
if 'table' in sdata.tables and 'metadata' not in sdata.tables:
  sdata['metadata'] = sdata.tables['table']
  del sdata.tables['table']

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
assert len(sdata_transcripts["transcripts"]["cell_id"].unique()) >= adata.n_obs, "Number of cells in transcripts must be at least the number of cells in the adata"

##################
# Aggregate data #
##################
sdata["transcripts"] = sdata_transcripts["transcripts"]
adata.obs['passed_QC'] = adata_qc_col.obs['passed_QC']
sdata['counts'] = adata

#######################
# Compute cell shapes #
#######################
print('Computing cell boundaries from transcripts using convex hulls', flush=True)
transcripts_df = sdata_transcripts["transcripts"].compute()
transcripts_assigned = transcripts_df[transcripts_df["cell_id"] != 0]
cell_shapes = transcripts_assigned.groupby("cell_id")[["x", "y"]].apply(
  lambda g: MultiPoint(list(zip(g["x"], g["y"]))).convex_hull
)
geo_df = gpd.GeoDataFrame(geometry=cell_shapes)
geo_df = sopa.shapes.to_valid_polygons(geo_df)
transformations = copy_transformations(sdata_transcripts["transcripts"])
sdata["cell_boundaries"] = ShapesModel.parse(geo_df, transformations=transformations)

#################
# Write output #
#################
sdata.write(par['output'])

