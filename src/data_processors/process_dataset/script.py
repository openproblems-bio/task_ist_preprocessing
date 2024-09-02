import scanpy as sc
import spatialdata as sd

### VIASH START
par = {
  "input_sc": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_sp": "resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
  "output_sc": "resources_test/preprocessing_imagingbased_st/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_sp": "resources_test/preprocessing_imagingbased_st/2023_10x_mouse_brain_xenium/dataset.zarr"
}
### VIASH END

# Load the single-cell data
input_sc = sc.read(par["input_sc"])

# Load the spatial data
input_sp = sd.read_zarr(par["input_sp"])

# Process if need be
output_sc = input_sc
output_sp = input_sp

# Save the single-cell data
output_sc.write_h5ad(par["output_sc"])

# Save the spatial data
output_sp.write(par["output_sp"])
