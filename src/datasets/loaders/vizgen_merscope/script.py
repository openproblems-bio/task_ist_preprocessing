# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

import os
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import h5py
import pandas as pd
import spatialdata as sd
import spatialdata_io as sdio
from shapely import MultiPolygon
from shapely.geometry import Polygon

## VIASH START
par = {
    "input": "path/to/HumanBreastCancerPatient1",
    "segmentation_id": ["cell"],
    "output": "output.zarr",
    "dataset_id": "vizgen_merscope/2022_vizgen_human_breast_cancer_merfish/rep1",
    "dataset_name": "value",
    "dataset_url": "https://vizgen.com/data-release-program/",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "human",
}


## VIASH END

assert ("cell" in par["segmentation_id"]) and (len(par["segmentation_id"]) == 1), (
    "Currently cell labels are definitely assigned in this script. And merscope does not provide other segmentations."
)

t0 = datetime.now()

print(datetime.now() - t0, "Starting vizgen_merscope preprocessing", flush=True)
print(
    datetime.now() - t0,
    "Parameters:",
    {k: par[k] for k in ["input", "output", "dataset_id", "dataset_organism"]},
    flush=True,
)


def read_boundary_hdf5(folder):
    print(datetime.now() - t0, "Convert boundary hdf5 to parquet", flush=True)
    all_boundaries = {}
    boundaries = None
    hdf5_files = os.listdir(folder + "/cell_boundaries/")
    n_files = len(hdf5_files)
    incr = n_files // 15
    for _, i in enumerate(hdf5_files):
        if (_ % incr) == 0:
            print(datetime.now() - t0, f"\tProcessed {_}/{n_files}", flush=True)
        with h5py.File(folder + "/cell_boundaries/" + i, "r") as f:
            for key in f["featuredata"].keys():
                if boundaries is not None:
                    boundaries.loc[key] = MultiPolygon(
                        [
                            Polygon(
                                f["featuredata"][key]["zIndex_3"]["p_0"]["coordinates"][
                                    ()
                                ][0]
                            )
                        ]
                    )  # doesn't matter which zIndex we use, MultiPolygon to work with read function in spatialdata-io
                else:
                    # boundaries = gpd.GeoDataFrame(index=[key], geometry=MultiPolygon([Polygon(f['featuredata'][key]['zIndex_3']['p_0']['coordinates'][()][0])]))
                    boundaries = gpd.GeoDataFrame(
                        geometry=gpd.GeoSeries(
                            MultiPolygon(
                                [
                                    Polygon(
                                        f["featuredata"][key]["zIndex_3"]["p_0"][
                                            "coordinates"
                                        ][()][0]
                                    )
                                ]
                            ),
                            index=[key],
                        )
                    )
            all_boundaries[i] = boundaries
            boundaries = None
    all_concat = pd.concat(list(all_boundaries.values()))
    all_concat = all_concat[
        ~all_concat.index.duplicated(keep="first")
    ]  # hdf5 can contain duplicates with same cell_id and position, removing those
    all_concat.rename_geometry(
        "geometry_renamed", inplace=True
    )  # renaming to make it compatible with spatialdata-io
    all_concat["EntityID"] = (
        all_concat.index
    )  # renaming to make it compatible with spatialdata-io
    all_concat["ZIndex"] = 0  # adding to make it compatible with spatialdata-io
    all_concat.to_parquet(folder + "/cell_boundaries.parquet")

    print(
        datetime.now() - t0,
        f"Wrote parquet: {folder}/cell_boundaries.parquet",
        flush=True,
    )

    count_path = f"{folder}/cell_by_gene.csv"
    obs_path = f"{folder}/cell_metadata.csv"

    data = pd.read_csv(count_path, index_col=0)
    obs = pd.read_csv(obs_path, index_col=0)

    data.index = obs.index.astype(str)  # data index in old software is range(n_obs)
    data.index.name = "cell"  # renaming to make it compatible with spatialdata-io
    obs.index.name = "EntityID"  # renaming to make it compatible with spatialdata-io
    data.to_csv(count_path)
    obs.to_csv(obs_path)

    print(
        datetime.now() - t0,
        f"Updated count/obs CSVs: {count_path}, {obs_path}",
        flush=True,
    )


RAW_DATA_DIR = Path(par["input"])

# If the cell polygons are in the old format (cell_boundaries/*.hdf5 instead of cell_boundaries.parquet) the raw data
# needs to be modified for the spatialdata-io loader
# (see: https://github.com/scverse/spatialdata-io/issues/71#issuecomment-1741995582)
if not (RAW_DATA_DIR / "cell_boundaries.parquet").exists():
    print(
        datetime.now() - t0, "Old boundary format detected, converting...", flush=True
    )
    read_boundary_hdf5(str(RAW_DATA_DIR))
    print(datetime.now() - t0, "Boundary conversion finished", flush=True)
else:
    print(
        datetime.now() - t0, "Boundary parquet found, skipping conversion", flush=True
    )


# Generate spatialdata.zarr

#########################################
# Convert raw files to spatialdata zarr #
#########################################
print(datetime.now() - t0, "Convert raw files to spatialdata zarr", flush=True)

slide_name = "slide"

print(datetime.now() - t0, "Calling spatialdata_io.merscope loader...", flush=True)
sdata = sdio.merscope(
    RAW_DATA_DIR,
    vpt_outputs=None,
    z_layers=3,
    region_name=None,
    slide_name=slide_name,
    backend=None,
    transcripts=True,
    cells_boundaries=True,
    cells_table=True,
    mosaic_images=True,
    # imread_kwargs=mappingproxy({}),
    # image_models_kwargs=mappingproxy({})
)
print(datetime.now() - t0, "Loader finished", flush=True)
try:
    print(
        datetime.now() - t0,
        "Loaded elements: " + ", ".join(list(sdata.keys())),
        flush=True,
    )
except Exception as e:
    print(
        datetime.now() - t0,
        f"Loaded spatialdata object (could not list keys): {type(e).__name__}: {e}",
        flush=True,
    )

###############
# Rename keys #
###############
print(datetime.now() - t0, "Rename keys", flush=True)

name = slide_name + "_" + RAW_DATA_DIR.name

elements_renaming_map = {
    f"{name}_z3": "morphology_mip",  # TODO: that is actually not the morphology_mip, i.e. either we should rename the label later, or we should actually project over z. But we also want to have 3d at some point anyway
    f"{name}_transcripts": "transcripts",
    f"{name}_polygons": "cell_boundaries",
    #"table": "metadata",
}

for old_key, new_key in elements_renaming_map.items():
    sdata[new_key] = sdata[old_key]
    del sdata[old_key]

print(datetime.now() - t0, "Renamed elements", flush=True)

# Rename transcript column
sdata["transcripts"] = sdata["transcripts"].rename(columns={"global_z": "z", "transcript_id": "ensembl_id"})#, "gene": "feature_name"})
if "gene" in sdata["transcripts"].columns: 
    # No idea why, but somehow dask dataframe renaming for the 'gene' column ends up in a key error when assigning it to sdata["transcripts"].
    # update: see https://github.com/scverse/spatialdata/issues/996
    sdata["transcripts"]["feature_name"] = sdata["transcripts"]["gene"]
    del sdata["transcripts"]["gene"]
    sdata['transcripts'].attrs["spatialdata_attrs"]["feature_key"] = "feature_name"
print(datetime.now() - t0, "Renamed transcripts column 'global_z' -> 'z' and 'gene' -> 'feature_name' and 'transcript_id' -> 'ensembl_id'", flush=True)

print(datetime.now() - t0, "Columns in sdata['transcripts']:", sdata["transcripts"].columns, flush=True)

#########################################
# Throw out all channels except of DAPI #
#########################################
print(datetime.now() - t0, "Throw out all channels except of DAPI", flush=True)

# TODO: in the future we want to keep PolyT and Cellbound1/2/3 stains. Note however, that somehow saving or plotting the sdata fails when
#       these channels aren't excluded, not sure why...
sdata["morphology_mip"] = sdata["morphology_mip"].sel(c=["DAPI"])
print(datetime.now() - t0, "Selected DAPI channel", flush=True)

#################################
# Get cell labels from polygons #
#################################
print(datetime.now() - t0, "Get cell labels from polygons", flush=True)

# TODO: Just note that currently the rasterize function has a bug, this error is small though with the given spatial resolution.
#       Check https://github.com/scverse/spatialdata/issues/165 for updates on this bug.
# NOTE: we need to iteratively rasterize (see here: https://github.com/scverse/spatialdata/issues/987)
img_extent = sd.get_extent(sdata["morphology_mip"])

N = 65535
n_cells = len(sdata["cell_boundaries"])
n_iter = n_cells // N + bool(n_cells % N)

rasterize_args = {
    "min_coordinate": [int(img_extent["x"][0]), int(img_extent["y"][0])],
    "max_coordinate": [int(img_extent["x"][1]), int(img_extent["y"][1])],
    "target_coordinate_system": "global",
    "target_unit_to_pixels": 1,
    "return_regions_as_labels": True,
}

for i in range(n_iter):
    print(
        datetime.now() - t0,
        f"Rasterizing iteration {i + 1}/{n_iter} (cells {i * N}..{min((i + 1) * N, n_cells)})",
        flush=True,
    )
    labels_image_ = sd.rasterize(
        sdata["cell_boundaries"].iloc[i * N : min((i + 1) * N, n_cells)],
        ["x", "y"],
        **rasterize_args,
    )
    if i == 0:
        labels_image = labels_image_
    else:
        labels_image.values[labels_image_.values > 0] = labels_image_.values[
            labels_image_.values > 0
        ]

print(datetime.now() - t0, "Rasterization finished", flush=True)
try:
    print(datetime.now() - t0, f"Label image shape: {labels_image.shape}", flush=True)
except Exception as e:
    print(datetime.now() - t0, f"Could not access label image shape: {e}", flush=True)

sdata["cell_labels"] = labels_image

del sdata["cell_labels"].attrs["label_index_to_category"]


##############################
# Add info to metadata table #
##############################
print(datetime.now() - t0, "Add metadata to table", flush=True)

# TODO: values as input variables
for key in [
    "dataset_id",
    "dataset_name",
    "dataset_url",
    "dataset_reference",
    "dataset_summary",
    "dataset_description",
    "dataset_organism",
    "segmentation_id",
]:
    sdata["table"].uns[key] = par[key]

print(datetime.now() - t0, "Metadata updated", flush=True)

#########
# Write #
#########
print(datetime.now() - t0, f"Writing to {par['output']}", flush=True)

sdata.write(par["output"])

print(datetime.now() - t0, f"Write completed: {par['output']}", flush=True)
print(datetime.now() - t0, "Done", flush=True)
