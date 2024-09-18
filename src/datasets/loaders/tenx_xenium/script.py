import spatialdata as sd
import anndata as ad
from spatialdata_io import xenium
import shutil
import os
import zipfile
import tempfile

## VIASH START
par = {
    "input": "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs",
    "replicate_id": "rep1",
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
    "output": "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zarr"
}
## VIASH END


# if input is a zip, extract it to a temporary folder
par_input = par["input"]
with tempfile.TemporaryDirectory() as tmpdirname:
    if zipfile.is_zipfile(par_input):
        print("Extracting input zip", flush=True)
        with zipfile.ZipFile(par_input, "r") as zip_ref:
            zip_ref.extractall(tmpdirname)
            par_input = tmpdirname

    # read the data
    sdata = xenium(
        path=par_input,
        n_jobs=meta["cpus"] or 1,
        cells_boundaries=True,
        nucleus_boundaries=True,
        morphology_focus=True,
        cells_as_circles=False,
    )

    # remove morphology_focus
    _ = sdata.images.pop("morphology_focus")

    print("Add uns to table", flush=True)
    new_uns = {
        "dataset_id": par["dataset_id"],
        "dataset_name": par["dataset_name"],
        "dataset_url": par["dataset_url"],
        "dataset_reference": par["dataset_reference"],
        "dataset_summary": par["dataset_summary"],
        "dataset_description": par["dataset_description"],
        "dataset_organism": par["dataset_organism"],
        "replicate_id": par["replicate_id"],
        "segmentation_id": par["segmentation_id"],
    }
    for key, value in new_uns.items():
        sdata.tables["table"].uns[key] = value

    print(f"Output: {sdata}", flush=True)

    print(f"Writing to '{par['output']}'", flush=True)
    if os.path.exists(par["output"]):
        shutil.rmtree(par["output"])

    sdata.write(par["output"])
