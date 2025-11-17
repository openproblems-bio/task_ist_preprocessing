

"""
Notes about data structure:

The directory structure that we need, looks as follows (expected by the sopa.io.cosmx function):
├── CellStatsDir ( = `DATA_DIR`)
│   ├── CellLabels
│   │    ├── CellLabels_F001.tif
│   │    ├── ...
│   │    └── CellLabels_F<last_fov_id>.tif
│   ├── Morphology2D
│   │    ├── <some_id>_F001.tif
│   │    ├── ...
│   │    └── <some_id>_F<last_fov_id>.tif
│   ├── <dataset_id>_exprMat_file.csv
│   ├── <dataset_id>_fov_positions_file.csv
│   ├── <dataset_id>_metadata_file.csv
│   ├── <dataset_id>_tx_file.csv
│   └── <dataset_id>-polygons.csv
├── (AnalysisResults)
└── (RunSummary)

The NSCLC samples download two files:
1. <sample>+SMI+Flat+data.tar.gz --> decompressed:

└── <sample>/<sample>-Flat_files_and_images ( = `DATA_DIR`)
    ├── CellLabels
    │    ├── CellLabels_F001.tif
    │    ├── ...
    │    └── CellLabels_F<last_fov_id>.tif
    ├── <dataset_id>_exprMat_file.csv
    ├── <dataset_id>_fov_positions_file.csv
    ├── <dataset_id>_metadata_file.csv
    ├── <dataset_id>_tx_file.csv
    └── NOTE: missing: <dataset_id>-polygons.csv
    
2. <sample>+RawMorphologyImages.tar.gz --> decompressed:

└── <sample> 2/<sample>-RawMorphologyImages ( = `MORPHOLOGY_DIR`)
    ├── <some_id>_F001_Z001.TIF
    ├── <some_id>_F001_Z002.TIF
    ├── ...
    ├── <some_id>_F001_Z008.TIF
    ├── ...
    └── <some_id>_F<last_fov_id>_Z008.TIF

The morphology images need to be moved to DATA_DIR / "Morphology2D". Also, we only take the 4th plane (Z004) of each image
and rename e.g. <some_id>_F001_Z004.TIF to <some_id>_F001.tif.


"""

import os
import shutil
import tarfile
from pathlib import Path
from datetime import datetime
import spatialdata as sd
import sopa

## VIASH START
par = {
    "sample": "Lung9_Rep1",
    "segmentation_id": ["cell"],
    "output": "output.zarr",
    "dataset_id": "bruker_cosmx/bruker_human_nsclc_cosmx/lung9_rep1",
    "dataset_name": "value",
    "dataset_url": "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "human",
}
meta = {
    #"temp_dir": "./temp/datasets/bruker_cosmx",
    "temp_dir": "/Volumes/T7/G3_temp/datasets/cosmx/test_folder"
}

## VIASH END

assert ("cell" in par["segmentation_id"]) and (len(par["segmentation_id"]) == 1), "Currently cell labels are definitely assigned in this script. And cosmx does not provide other segmentations."

t0 = datetime.now()

# Download urls for sample
URL = f"https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/{par['sample']}/{par['sample']}+SMI+Flat+data.tar.gz"
URL_MORPHOLOGY = f"https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/{par['sample']}/{par['sample']}+RawMorphologyImages.tar.gz"

# Define temp dir and file names
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_NAME = TMP_DIR / (par["sample"] + "+SMI+Flat+data.tar.gz")
FILE_NAME_MORPHOLOGY = TMP_DIR / (par["sample"] + "+RawMorphologyImages.tar.gz")
DATA_DIR = TMP_DIR / par["sample"] / (par["sample"] + "-Flat_files_and_images")
MORPHOLOGY_DIR = TMP_DIR / par["sample"] / (par["sample"] + "-RawMorphologyImages")

# Download flat files and cell labels
print(datetime.now() - t0, "Download flat files and cell labels", flush=True)
os.system(f"wget {URL} -O '{FILE_NAME}'")

# Extract tar.gz files
print(datetime.now() - t0, "Extract tar.gz of flat files and cell labels", flush=True)
with tarfile.open(FILE_NAME, "r:gz") as tar:
    tar.extractall(TMP_DIR)

# Check that flat files are present
FLAT_FILES_ENDINGS = ["_exprMat_file.csv", "_fov_positions_file.csv", "_metadata_file.csv", "_tx_file.csv"] #, "polygons.csv"]
for ending in FLAT_FILES_ENDINGS:
    if any(f.endswith(ending) for f in os.listdir(DATA_DIR)):
        print(f"Flat file with ending {ending} is present", flush=True)
    else:
        print(f"Flat file with ending {ending} is missing", flush=True)

# Download image files
print(datetime.now() - t0, "Download image files", flush=True)
os.system(f"wget {URL_MORPHOLOGY} -O '{FILE_NAME_MORPHOLOGY}'")

# Extract tar.gz files
print(datetime.now() - t0, "Extract tar.gz of image files", flush=True)
with tarfile.open(FILE_NAME_MORPHOLOGY, "r:gz") as tar:
    tar.extractall(TMP_DIR)

# Move image files to DATA_DIR / "Morphology2D"
print(datetime.now() - t0, f"Move image files of plane Z004 to {DATA_DIR / 'Morphology2D'} and rename.", flush=True)
(DATA_DIR / "Morphology2D").mkdir(parents=True, exist_ok=True)
for image_file in MORPHOLOGY_DIR.glob("*_Z004.TIF"):
    shutil.move(image_file, DATA_DIR / "Morphology2D" / image_file.name.replace("_Z004.TIF", ".tif"))



#########################################
# Convert raw files to spatialdata zarr #
#########################################

# from pathlib import Path
# import sopa
# DATA_DIR = Path("/Volumes/T7/G3_temp/datasets/cosmx/test_folder/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images")

print(datetime.now() - t0, "Convert raw files to spatialdata zarr", flush=True)

# The tif images of the NSCLC samples miss metadata information, so we hardcode the channels here.
def fixed_get_morphology_coords(images_dir: Path) -> list[str]:
    return ["other1", "other2", "other3", "other4", "DAPI"]

sopa.io.reader.cosmx._get_morphology_coords = fixed_get_morphology_coords

sdata = sopa.io.cosmx(
    DATA_DIR, 
    dataset_id=None, 
    fov=None, 
    read_proteins=False, 
    cells_labels=True, 
    cells_table=True, 
    cells_polygons=False, 
    flip_image=False
)

# Retrieve polygons from segmentation (no polygons file in flat files present for NSCLC samples)
sdata["cells_polygons"] = sd.to_polygons(sdata["stitched_labels"])


###############
# Rename keys #
###############
print(datetime.now() - t0, "Rename keys", flush=True)

elements_renaming_map = {
    "stitched_image"     : "morphology_mip", 
    "stitched_labels"    : "cell_labels",
    "points"             : "transcripts",
    "cells_polygons"     : "cell_boundaries",
    #"table"              : "metadata",
}

for old_key, new_key in elements_renaming_map.items():
    sdata[new_key] = sdata[old_key]
    del sdata[old_key]

# Rename transcript column (somehow overwriting the 'target' column leads to an error, so instead we add a duplicate with the right name)
#sdata['transcripts'] = sdata['transcripts'].rename(columns={"global_cell_id":"cell_id", "target":"feature_name"})
sdata['transcripts'] = sdata['transcripts'].rename(columns={"global_cell_id":"cell_id"})
sdata['transcripts']["feature_name"] = sdata['transcripts']["target"]

#########################################
# Throw out all channels except of DAPI #
#########################################
print(datetime.now() - t0, "Throw out all channels except of 'DAPI'", flush=True)

# TODO: in the future we want to keep PolyT and Cellbound1/2/3 stains. Note however, that somehow saving or plotting the sdata fails when
#       these channels aren't excluded, not sure why...
sdata["morphology_mip"] = sdata["morphology_mip"].sel(c=["DAPI"])


##############################
# Add info to metadata table #
##############################
print(datetime.now() - t0, "Add metadata to table", flush=True)

#TODO: values as input variables
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism", "segmentation_id"]:
    sdata["table"].uns[key] = par[key]

#########
# Write #
#########
print(datetime.now() - t0, f"Writing to {par['output']}", flush=True)

sdata.write(par["output"])

print(datetime.now() - t0, "Done", flush=True)