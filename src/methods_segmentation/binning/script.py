import txsim as tx
import numpy as np
import os
import yaml
import spatialdata as sd
import anndata as ad
import shutil
import numpy as np
from spatialdata.models import Labels2DModel
import xarray as xr



def convert_to_lower_dtype(arr):
    max_val = arr.max()
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64

    return arr.astype(new_dtype)

## VIASH START
par = {
  "input": "../task_ist_preprocessing/resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
  "output": "segmentation.zarr"
}

## VIASH END

hyperparameters = par.copy()

hyperparameters = {k:(v if v != "None" else None) for k,v in hyperparameters.items()}
del hyperparameters['input']
del hyperparameters['output']

sdata = sd.read_zarr(par["input"])
image = sdata['morphology_mip']['scale0'].image.compute().to_numpy()
transformation = sdata['morphology_mip']['scale0'].image.transform.copy()

sd_output = sd.SpatialData()
image = sdata['morphology_mip']['scale0'].image.compute().to_numpy()
transformation = sdata['morphology_mip']['scale0'].image.transform.copy()
img_arr = tx.preprocessing.segment_binning(image[0], hyperparameters['bin_size'])   ### TOdo find the optimal bin_size
image = convert_to_lower_dtype(img_arr)
data_array = xr.DataArray(image, name=f'segmentation', dims=('y', 'x'))
parsed_data = Labels2DModel.parse(data_array, transformations=transformation)
sd_output.labels['segmentation'] = parsed_data

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sd_output.write(par["output"])

