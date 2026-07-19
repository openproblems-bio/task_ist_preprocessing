import spatialdata as sd
import anndata as ad
import os
import shutil
import pandas as pd

## VIASH START
par = {
  "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
  "labels_key": "cell_labels",
  "output": "segmentation.zarr",
}
meta = {
  "name": "segmentation"
}
## VIASH END

print("Reading input files", flush=True)
sdata = sd.read_zarr(par["input"])

assert par["labels_key"] in sdata.labels, f"Key '{par['labels_key']}' not found in input data."

print(f"Copy segmentation from '{par['labels_key']}'", flush=True)
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
sdata_segmentation_only = sd.SpatialData(
  labels={
    "segmentation": sdata[par["labels_key"]]
  },
  tables={
    "table": ad.AnnData(
      obs=obs,
      var=metadata.var[[]]
    )
  }
)

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_segmentation_only.write(par["output"])
