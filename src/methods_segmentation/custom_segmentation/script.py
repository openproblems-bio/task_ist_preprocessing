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
# Select only the columns that exist — Xenium provides cell_id and region,
# Vizgen uses different column names (or an empty obs) so we take what's available.
obs_cols = [c for c in ["cell_id", "region"] if c in metadata.obs.columns]
sdata_segmentation_only = sd.SpatialData(
  labels={
    "segmentation": sdata[par["labels_key"]]
  },
  tables={
    "table": ad.AnnData(
      obs=metadata.obs[obs_cols],
      var=metadata.var[[]]
    )
  }
)

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_segmentation_only.write(par["output"])
