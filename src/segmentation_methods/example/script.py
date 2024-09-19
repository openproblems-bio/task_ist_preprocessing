import spatialdata as sd

### VIASH START
par = {
  "input": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset_segmented.zarr"
}
### VIASH END

# Load the spatial data
input = sd.read_zarr(par["input"])

# Do something with input
output = input

# Save the spatial data
output.write(par["output"])
