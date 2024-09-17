import spatialdata as sd
import anndata as ad
from spatialdata_io import xenium
import shutil
import os
import zipfile
import tempfile

## VIASH START
par = {
    "input": [
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs",
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs",
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    ],
    "replicate_id": [
        "rep1",
        "rep2",
        "rep3",
    ],
    "segmentation_id": [
        "cell",
        "nucleus",
    ],
    "output": "output.zarr",
    "dataset_id": "value",
    "dataset_name": "value",
    "dataset_url": "value",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "value",
}
## VIASH END

assert len(par["input"]) == len(par["replicate_id"]), "Length of 'input' and 'replicate_id' must be the same."

sdatas = []

for i, input in enumerate(par["input"]):
    replicate_id = par["replicate_id"][i]
    print(f"Processing replicate '{replicate_id}'", flush=True)

    # if input is a zip, extract it to a temporary folder
    with tempfile.TemporaryDirectory() as tmpdirname:
        if zipfile.is_zipfile(input):
            print("Extracting input zip", flush=True)
            with zipfile.ZipFile(input, "r") as zip_ref:
                zip_ref.extractall(tmpdirname)
                input = tmpdirname

        # read the data
        sdata = xenium(
            path=input,
            n_jobs=8,
            cells_boundaries=True,
            nucleus_boundaries=True,
            morphology_focus=True,
            cells_as_circles=False,
        )

    # rename coordinate system
    sdata.rename_coordinate_systems({"global": replicate_id + "_global"})
    
    # rename images
    sdata.images[replicate_id + "_image"] = sdata.images.pop("morphology_mip")

    # remove morphology_focus
    _ = sdata.images.pop("morphology_focus")

    # rename labels
    sdata.labels[replicate_id + "_cell"] = sdata.labels.pop("cell_labels")
    sdata.labels[replicate_id + "_nucleus"] = sdata.labels.pop("nucleus_labels")

    # rename points
    sdata.points[replicate_id + "_transcripts"] = sdata.points.pop("transcripts")

    # rename shapes
    sdata.shapes[replicate_id + "_cell_boundaries"] = sdata.shapes.pop("cell_boundaries")
    sdata.shapes[replicate_id + "_nucleus_boundaries"] = sdata.shapes.pop("nucleus_boundaries")

    # rename tables
    sdata.tables[replicate_id + "_cell_table"] = sdata.tables.pop("table")

    sdatas.append(sdata)

print("Concatenate sdatas", flush=True)
sdata = sd.concatenate(sdatas)

print("Add metadata table", flush=True)
sdata.tables["metadata"] = ad.AnnData(
    uns={
        "dataset_id": par["dataset_id"],
        "dataset_name": par["dataset_name"],
        "dataset_url": par["dataset_url"],
        "dataset_reference": par["dataset_reference"],
        "dataset_summary": par["dataset_summary"],
        "dataset_description": par["dataset_description"],
        "dataset_organism": par["dataset_organism"],
        "variables": {
            "replicate_id": par["replicate_id"],
            "segmentation_id": par["segmentation_id"],
        }
    }
)

print(f"Output: {sdata}", flush=True)

print(f"Writing to '{par['output']}'", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])

sdata.write(par["output"])
