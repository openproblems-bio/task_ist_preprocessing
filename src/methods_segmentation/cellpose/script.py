import txsim as tx
#from pathlib import Path
#import tifffile
#import squidpy as sq
#from scipy import ndimage
#import skimage.io
#import skimage.measure
#import skimage.segmentation
import numpy as np
import txsim as tx 
#import argparse
import os
import yaml
import spatialdata as sd
import anndata as ad
import shutil
import numpy as np
from spatialdata.models import Labels2DModel
import xarray as xr
import datatree as dt


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
  "output": "segmentation.zarr",
  'hyperparams': [] ##############################how to take from viash
}

## VIASH END

sdata = sd.read_zarr(par["input"])
image = sdata['rep1_image']['scale0'].image.compute().to_numpy()
transformation = image.transform.copy()
transformation['global'] = transformation.pop('rep1_global')
image = convert_to_lower_dtype(image)

sd_output = sd.SpatialData()
scales = [2000, 1000, 500, 250, 125]
downsampled_arrays = {}
for idx, scale in enumerate(scales):
    img_arr = tx.preprocessing.segment_cellpose(image, hyperparams)  
    data_array = xr.DataArray(img_arr, name=f'segmentation_scale{scale}', dims=('y', 'x'))
    parsed_data = Labels2DModel.parse(data_array, transformations=transformation)
    downsampled_arrays[f'scale{idx}'] = parsed_data
tree = dt.DataTree()
for scale_key, parsed_data in downsampled_arrays.items():
    tree[scale_key] = dt.DataTree(parsed_data)

sd_output.labels['segmentation'] = tree


print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sd_output.write(par["output"])

