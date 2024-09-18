import spatialdata as sd
import shutil
import os

## VIASH START
par = {
    "input": "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zarr",
    "output": "resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
    "min_x": 10000,
    "max_x": 12000,
    "min_y": 10000,
    "max_y": 12000,
}
## VIASH END

print("Reading input", flush=True)
sdata = sd.read_zarr(par["input"])

print(f"  Cropping to region: {par['min_y']}:{par['max_y']} x {par['min_x']}:{par['max_x']}", flush=True)
sdata_output = sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[par["min_x"], par["min_y"]],
    max_coordinate=[par["max_x"], par["max_y"]],
    target_coordinate_system="global",
    filter_table=True,
)

print(f"Output: {sdata_output}", flush=True)

# write the output
print(f"Writing output to '{par['output']}'", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sdata_output.write(par["output"])
