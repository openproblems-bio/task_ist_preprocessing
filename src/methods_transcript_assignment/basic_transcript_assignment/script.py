import numpy as np
import xarray as xr
import dask
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
transcripts = sd.transform(sdata[par['transcripts_key']], to_coordinate_system=par['coordinate_system'])

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
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        label_image[y_coords, x_coords], chunks=tuple(sdata[par['transcripts_key']].map_partitions(len).compute())
    ), 
    index=sdata[par['transcripts_key']].index
)
sdata[par['transcripts_key']]["cell_id"] = cell_id_dask_series 

#create new .obs for cells based on the segmentation output (corresponding with the transcripts 'cell_id')
unique_cells = np.unique(cell_id_dask_series)

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
    "transcripts": sdata[par['transcripts_key']]
  },
  tables={
    "table": ad.AnnData(
      obs=pd.DataFrame(cell_id_col),
      var=sdata.tables["table"].var[[]]
    )
  }
)

print('Write transcripts with cell ids', flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_transcripts_only.write(par['output'])
