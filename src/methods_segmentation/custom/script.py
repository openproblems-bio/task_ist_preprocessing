import spatialdata as sd
import anndata as ad
import os
import shutil

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
sdata_segmentation_only = sd.SpatialData(
  labels={
    "segmentation": sdata[par["labels_key"]]
  },
  tables={
    "table": ad.AnnData(
      obs=sdata.tables["table"].obs[["cell_id", "region"]],
      var=sdata.tables["table"].var[[]]
    )
  }
)

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_segmentation_only.write(par["output"])
