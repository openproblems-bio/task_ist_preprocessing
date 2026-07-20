import numpy as np
import xarray as xr
import dask
import dask.dataframe as dd
import spatialdata as sd
import anndata as ad
import pandas as pd
import os
import shutil

## VIASH START
par = {
  'input_ist': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
  'input_segmentation': 'resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr',
  'transcripts_key': 'transcripts',
  'coordinate_system': 'global',
  'output': 'basic_assigned_transcripts.zarr',
}
meta = {
  'name': 'basic'
}
## VIASH END

assert par['input_segmentation'], 'Segmentation input is required for this assignment method.'

print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])

# Check if coordinate system is available in input data
transcripts_coord_systems = sd.transformations.get_transformation(sdata[par["transcripts_key"]], get_all=True).keys()
assert par['coordinate_system'] in transcripts_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."
segmentation_coord_systems = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).keys()
assert par['coordinate_system'] in segmentation_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."

print('Transforming transcripts coordinates', flush=True)
# Multi-partition parquet files each start with a 0-based index, producing duplicate index
# values in the combined dask DataFrame. sd.transform() internally creates a pd.Series with
# index=transformed.index; when that dask index is computed it triggers an assign expression
# that fails on duplicate/lazy indices. Fix: materialize to pandas and rebuild as a single
# dask partition with a clean RangeIndex before transforming.
# The original sdata[transcripts_key] is left unchanged so lines below remain consistent.
transcripts_input = sdata[par['transcripts_key']]
transcripts_reset = dd.from_pandas(transcripts_input.compute().reset_index(drop=True), npartitions=1)
transcripts_reset.attrs.update(transcripts_input.attrs)
transcripts = sd.transform(transcripts_reset, to_coordinate_system=par['coordinate_system'])

# In case of a translation transformation of the segmentation (e.g. crop of the data), we need to adjust the transcript coordinates
trans = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True)[par['coordinate_system']].inverse()
transcripts = sd.transform(transcripts, trans, par['coordinate_system'])

print('Assigning transcripts to cell ids', flush=True)
y_coords = transcripts.y.compute().to_numpy(dtype=np.int64)
x_coords = transcripts.x.compute().to_numpy(dtype=np.int64)
if isinstance(sdata_segm["segmentation"], xr.DataTree):
    label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy() 
else:
     label_image = sdata_segm["segmentation"].to_numpy()
# Assign cell ids directly to transcripts_reset (clean-index single-partition dask DataFrame).
# Using sdata[par['transcripts_key']] here would reintroduce the duplicate parquet index,
# causing the same "cannot reindex on an axis with duplicate labels" error at write time.
# Clamp coords to the label-image bounds: transcripts at the crop boundary can
# round a few pixels past the raster edge (see get_crop_coords in
# process_dataset). Matches segger's handling; edge/background at the border.
y_coords = np.clip(y_coords, 0, label_image.shape[0] - 1)
x_coords = np.clip(x_coords, 0, label_image.shape[1] - 1)
cell_ids = label_image[y_coords, x_coords]
transcripts_reset["cell_id"] = pd.Series(cell_ids)

#create new .obs for cells based on the segmentation output (corresponding with the transcripts 'cell_id')
unique_cells = np.unique(cell_ids)

# check if a '0' (noise/background) cell is in cell_id and remove
zero_idx = np.where(unique_cells == 0)
if len(zero_idx[0]): unique_cells=np.delete(unique_cells, zero_idx[0][0])

#transform into pandas series and check
cell_id_col = pd.Series(unique_cells, name='cell_id', index=unique_cells)
assert 0 not in cell_id_col, "Found '0' in cell_id column of assingment output cell matrix"

# TODO: Also take care of the following cases:
# - segmentation 3D, transcripts 3D
# - segmentation 3D, transcripts 2D
# - segmentation 2D, transcripts 3D

print('Subsetting to transcripts cell id data', flush=True)
sdata_transcripts_only = sd.SpatialData(
  points={
    "transcripts": transcripts_reset
  },
  tables={
    "table": ad.AnnData(
      obs=pd.DataFrame(cell_id_col),
    )
  }
)

print('Write transcripts with cell ids', flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_transcripts_only.write(par['output'])
