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

sd_output = sd.SpatialData()
image = sdata['image']['scale0'].image.compute().to_numpy()
transformation = sdata['image']['scale0'].image.transform.copy()
img_arr = tx.preprocessing.segment_watershed(image[0], hyperparameters) 
image = convert_to_lower_dtype(img_arr)
data_array = xr.DataArray(image, name=f'segmentation', dims=('y', 'x'))
parsed_data = Labels2DModel.parse(data_array, transformations=transformation)
sd_output.labels['segmentation'] = parsed_data

metadata = sdata.tables["metadata"]
# cell_id is required downstream. Standard Xenium exports carry an explicit
# "cell_id" column, but some exports (e.g. the Xenium WTA preview used for the
# Atera dataset) don't — there the per-cell identifier lives in the table's
# instance_key column (falling back to the obs index).
instance_key = metadata.uns.get("spatialdata_attrs", {}).get("instance_key")
if "cell_id" in metadata.obs.columns:
    cell_id = metadata.obs["cell_id"].values
elif instance_key and instance_key in metadata.obs.columns:
    cell_id = metadata.obs[instance_key].values
else:
    cell_id = metadata.obs.index.values
obs = metadata.obs[[]].copy()
obs["cell_id"] = cell_id
if "region" in metadata.obs.columns:
    obs["region"] = metadata.obs["region"].values
sd_output.tables['table'] = ad.AnnData(
      obs=obs,
      var=metadata.var[[]]
    )

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sd_output.write(par["output"])

