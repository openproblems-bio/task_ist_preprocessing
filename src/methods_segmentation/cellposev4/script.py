import anndata as ad
import numpy as np
import os
import shutil
import spatialdata as sd
import xarray as xr
from cellpose.models import CellposeModel
from spatialdata.models import Labels2DModel
import torch

# Check whether a GPU is available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('Using device:', device, flush=True)

## VIASH START
par = {
  "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
  "output": "segmentation.zarr",
  "diameter": 30.0,
  "flow_threshold": 0.0,
  "niter": 10,
  "min_size": -1,
  "resample": False,
}
meta = {
  "name": "cellposev4",
}
## VIASH END


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


print('Reading input', flush=True)
sdata = sd.read_zarr(par["input"])
image = sdata['image']['scale0'].image.compute().to_numpy()
transformation = sdata['image']['scale0'].image.transform.copy()

print('Initializing Cellpose model', flush=True)
model = CellposeModel(gpu=torch.cuda.is_available())

eval_params = {k: par[k] for k in ("diameter", "flow_threshold", "niter", "min_size", "resample") if par.get(k) is not None}
print('Running Cellpose segmentation with parameters:', eval_params, flush=True)
masks, _, _ = model.eval(image[0], progress=True, **eval_params)

print('Cellpose segmentation finished, post-processing results', flush=True)
masks = convert_to_lower_dtype(masks)

data_array = xr.DataArray(masks, name='segmentation', dims=('y', 'x'))
parsed_data = Labels2DModel.parse(data_array, transformations=transformation)

sd_output = sd.SpatialData()
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
