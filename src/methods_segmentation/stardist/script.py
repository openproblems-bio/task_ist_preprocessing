import os
import shutil
from pathlib import Path
import numpy as np
import xarray as xr
import spatialdata as sd
from csbdeep.utils import normalize
from stardist.models import StarDist2D



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
  "input": "./resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output": "./temp/stardist/segmentation.zarr",
  "model": "2D_versatile_fluo"
}

## VIASH END


# Read image and its transformation
sdata = sd.read_zarr(par["input"])
image = sdata['morphology_mip']['scale0'].image.compute().to_numpy()
transformation = sdata['morphology_mip']['scale0'].image.transform.copy()

# Segment image
# Load pretrained model
model = StarDist2D.from_pretrained(par['model'])
# Segment on normalized image 
labels, _ = model.predict_instances(normalize(image)[0,:,:]) # scale = None, **hyperparams)


# Create output
sd_output = sd.SpatialData()
labels = convert_to_lower_dtype(labels)
labels_array = xr.DataArray(labels, name=f'segmentation', dims=('y', 'x'))
parsed_labels = sd.models.Labels2DModel.parse(labels_array, transformations=transformation)
sd_output.labels['segmentation'] = parsed_labels

print("Writing output", flush=True)
Path(par["output"]).parent.mkdir(parents=True, exist_ok=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sd_output.write(par["output"])

