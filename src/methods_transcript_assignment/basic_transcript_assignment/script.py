import numpy as np
import dask
import spatialdata as sd
import anndata as ad
import os
import shutil

## VIASH START
par = {
  'input_ist': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
  'input_segmentation': 'resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr',
  'transcripts_key': 'transcripts',
  'coordinate_system': 'global',
  'output': 'assigned_transcripts.zarr',
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
label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy() #TODO: mabye this line needs generalization (DataTree vs DataArray)
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        label_image[y_coords, x_coords], chunks=tuple(sdata[par['transcripts_key']].map_partitions(len).compute())
    ), 
    index=sdata[par['transcripts_key']].index
)
sdata[par['transcripts_key']]["cell_id"] = cell_id_dask_series 

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
      obs=sdata.tables["table"].obs[["cell_id", "region"]],
      var=sdata.tables["table"].var[[]]
    )
  }
)

print('Write transcripts with cell ids', flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_transcripts_only.write(par['output'])
