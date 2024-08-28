import spatialdata as sd

## VIASH START
par = {
  "input": "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/full_dataset.zarr",
  "output": "resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
  "replicate_id": ["rep1", "rep2", "rep3"],
  "min_x": [10000, 10000, 10000],
  "max_x": [12000, 12000, 12000],
  "min_y": [10000, 10000, 10000],
  "max_y": [12000, 12000, 12000],
}
## VIASH END

sdata = sd.read_zarr(par["input"])

sdata_out = []

for i, replicate_id in enumerate(par["replicate_id"]):
  min_x = par["min_x"][i]
  max_x = par["max_x"][i]
  min_y = par["min_y"][i]
  max_y = par["max_y"][i]
  sdata_query = sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[min_x, min_y],
    max_coordinate=[max_x, max_y],
    target_coordinate_system=f"{replicate_id}_global",
  )
  sdata_out.append(sdata_query)

sdata_output = sd.concatenate(sdata_out)

sdata_output.write(par["output"])
