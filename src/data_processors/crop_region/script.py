import spatialdata as sd
import shutil
import os

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

print("Reading input", flush=True)
sdata = sd.read_zarr(par["input"])

print("Set aside metadata", flush=True)
sdata_metadata = sdata.tables.pop("metadata")

print("Extracting variables", flush=True)
replicate_ids = sdata_metadata.uns["variables"]["replicate_id"]
segmentation_ids = sdata_metadata.uns["variables"]["segmentation_id"]

# process the replicates
sdata_out = []

for i, replicate_id in enumerate(par["replicate_id"]):
    print(f"Processing replicate '{replicate_id}'", flush=True)

    min_x = par["min_x"][i]
    max_x = par["max_x"][i]
    min_y = par["min_y"][i]
    max_y = par["max_y"][i]

    print(f"  Cropping to region: {min_x}, {min_y}, {max_x}, {max_y}", flush=True)
    sdata_query = sdata.query.bounding_box(
        axes=["x", "y"],
        min_coordinate=[min_x, min_y],
        max_coordinate=[max_x, max_y],
        target_coordinate_system=f"{replicate_id}_global",
        filter_table=True,
    )

    # process the segmentations
    for segmentation_id in segmentation_ids:
        print(
            f"  Processing replicate '{replicate_id}' segmentation '{segmentation_id}'",
            flush=True,
        )
        shape_name = f"{replicate_id}_{segmentation_id}_boundaries"
        table_name = f"{replicate_id}_{segmentation_id}_table"

        if shape_name not in sdata_query.shapes:
            print(f"    Shape '{shape_name}' not found in sdata, skipping", flush=True)
            continue

        if table_name not in sdata.tables:
            print(f"    Table '{table_name}' not found in sdata, skipping", flush=True)
            continue

        obs_index = sdata_query.shapes[shape_name].index
        table_query = sdata.tables[table_name][obs_index]

        # add filtered table to sdata_query
        sdata_query.tables[table_name] = table_query

    # add filtered shapes to sdata_query
    sdata_out.append(sdata_query)

# concatenate the sdata objects
print("Concatenating sdata objects", flush=True)
sdata_output = sd.concatenate(sdata_out)

# add metadata back
print("Adding metadata back", flush=True)
sdata_output.tables["metadata"] = sdata_metadata.copy()

print(f"Output: {sdata_output}", flush=True)

# write the output
print(f"Writing output to '{par['output']}'", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sdata_output.write(par["output"])
