import spatialdata as sd

## VIASH START
par = {
  "input": "resources/datasets/10x_xenium/10x_fresh_frozen_mouse_brain_replicates/dataset.zarr",
  "output": "output.zarr",
  "replicate": ["rep1"],
  "min_x": [10000],
  "max_x": [12000],
  "min_y": [10000],
  "max_y": [12000]
}
## VIASH END

sdata = sd.read_zarr(par["input"])

for i, replicate_id in enumerate(par["replicate"]):
  min_x = par["min_x"][i]
  max_x = par["max_x"][i]
  min_y = par["min_y"][i]
  max_y = par["max_y"][i]
  sdata = sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[min_x, min_y],
    max_coordinate=[max_x, max_y],
    target_coordinate_system=f"{replicate_id}_global",
  )

sdata.write_zarr(par["output"])
