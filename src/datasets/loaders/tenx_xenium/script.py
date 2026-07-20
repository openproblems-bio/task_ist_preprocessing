import spatialdata as sd
import anndata as ad
from spatialdata_io import xenium
import shutil
import os
import zipfile
import tempfile

## VIASH START
par = {
    "input": "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hLiver_cancer_section_FFPE/Xenium_V1_hLiver_cancer_section_FFPE_outs.zip",
    "segmentation_id": [
        "cell",
        "nucleus",
    ],
    "dataset_id": "value",
    "dataset_name": "value",
    "dataset_url": "value",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "value",
    "output": "temp/datasets/10x_xenium/liver/liver.zarr"
}
meta = {
    "cpus": 1,
}

## VIASH END

# Download the data if it's a download url, extract the data if it's a zip file
par_input = par["input"]
with tempfile.TemporaryDirectory() as tmpdirname:
    if par_input.startswith("http"):
        print(f"Downloading data to {tmpdirname}", flush=True)
        file_name = par_input.split("/")[-1]
        # download the data
        os.system(f"wget {par['input']} -O {tmpdirname}/{file_name}")
        par_input = tmpdirname + "/" + file_name

    if zipfile.is_zipfile(par_input):
        print(f"Extracting input zip to {tmpdirname}", flush=True)
        with zipfile.ZipFile(par_input, "r") as zip_ref:
            zip_ref.extractall(tmpdirname)
        # find the directory containing the Xenium output files (may be nested)
        par_input = tmpdirname
        for root, dirs, files in os.walk(tmpdirname):
            if "cell_feature_matrix.h5" in files:
                par_input = root
                break

    # read the data
    sdata = xenium(
        path=par_input,
        n_jobs=meta["cpus"] or 1,
        cells_boundaries=True,
        nucleus_boundaries=True,
        morphology_focus=True,
        cells_as_circles=False,
    )

    # Keep exactly one morphology raster named "morphology_mip"; process_dataset
    # renames it to the API-required "image" element. Older Xenium ship both a
    # single-channel "morphology_mip" and a multi-channel "morphology_focus";
    # newer Xenium (Onboard Analysis v2+) ship only "morphology_focus" (channel 0
    # is DAPI, which the segmentation methods use via image[0]). Prefer the MIP,
    # otherwise fall back to morphology_focus so the dataset always has an image.
    if "morphology_mip" not in sdata.images and "morphology_focus" in sdata.images:
        sdata["morphology_mip"] = sdata["morphology_focus"]
    if "morphology_focus" in sdata.images:
        del sdata.images["morphology_focus"]

    print("Add uns to table", flush=True)
    new_uns = {
        "dataset_id": par["dataset_id"],
        "dataset_name": par["dataset_name"],
        "dataset_url": par["dataset_url"],
        "dataset_reference": par["dataset_reference"],
        "dataset_summary": par["dataset_summary"],
        "dataset_description": par["dataset_description"],
        "dataset_organism": par["dataset_organism"],
        "segmentation_id": par["segmentation_id"],
    }
    for key, value in new_uns.items():
        sdata.tables["table"].uns[key] = value

    print(f"Output: {sdata}", flush=True)

    print(f"Writing to '{par['output']}'", flush=True)
    if os.path.exists(par["output"]):
        shutil.rmtree(par["output"])

    print(f"Output: {sdata}")

    sdata.write(par["output"])
