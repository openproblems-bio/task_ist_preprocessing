import spatialdata as sd

### VIASH START
par = {
  "input": "resources_test/preprocessing_imagingbased_st/2023_10x_mouse_brain_xenium/dataset.zarr",
  "output": "resources_test/preprocessing_imagingbased_st/2023_10x_mouse_brain_xenium/dataset_segmented.zarr"
}
### VIASH END

# Load the spatial data
input = sd.read_zarr(par["input"])

# Do something with input
output = input

# Save the spatial data
output.write(par["output"])
