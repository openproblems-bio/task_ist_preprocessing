# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

import os
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import dask.array as da
import geopandas as gpd
import h5py
import numpy as np
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


def _element_names(sd_obj):
    """List element names across spatialdata versions.

    Newer spatialdata removed ``SpatialData.keys()``, so ``list(sdata.keys())`` raises
    ``AttributeError``. The per-type element dicts (``images``/``labels``/``points``/``shapes``/
    ``tables``) have been stable for far longer, so we collect names from those instead.
    """
    names = []
    for attr in ("images", "labels", "points", "shapes", "tables"):
        try:
            names.extend(getattr(sd_obj, attr).keys())
        except Exception:
            pass
    return names


def _cell_polygon(cell_grp):
    """Boundary for one cell: prefer zIndex_3, else any available z-plane.

    MERSCOPE stores a cell's boundary only at the z-planes where it appears, so keying on a
    single hardcoded z-index (the old ``["zIndex_3"]``) drops -- or KeyErrors on -- every cell
    not present at that plane. Returns None if the cell has no usable boundary.
    """
    z_keys = list(cell_grp.keys())
    z = "zIndex_3" if "zIndex_3" in z_keys else (z_keys[0] if z_keys else None)
    if z is None or "p_0" not in cell_grp[z] or "coordinates" not in cell_grp[z]["p_0"]:
        return None
    return MultiPolygon([Polygon(cell_grp[z]["p_0"]["coordinates"][()][0])])


def read_boundary_hdf5(folder):
    print(datetime.now() - t0, "Convert boundary hdf5 to parquet", flush=True)
    hdf5_files = [
        f for f in os.listdir(folder + "/cell_boundaries/") if f.endswith(".hdf5")
    ]
    n_files = len(hdf5_files)
    incr = max(1, n_files // 15)  # guard against ZeroDivision when <15 files
    geoms, index, n_skipped = [], [], 0
    for n, fname in enumerate(hdf5_files):
        if n % incr == 0:
            print(datetime.now() - t0, f"\tProcessed {n}/{n_files}", flush=True)
        with h5py.File(folder + "/cell_boundaries/" + fname, "r") as f:
            featuredata = f["featuredata"]
            for key in featuredata.keys():
                poly = _cell_polygon(featuredata[key])
                if poly is None:
                    n_skipped += 1
                    continue
                geoms.append(poly)
                index.append(key)
    # n_skipped is logged (never silently dropped) so a re-run reveals whether low boundary
    # coverage is a reader problem or a genuine property of the raw hdf5.
    print(
        datetime.now() - t0,
        f"Read {len(geoms)} cell boundaries from {n_files} files "
        f"({n_skipped} cells had no usable z-plane boundary)",
        flush=True,
    )
    # Build the GeoDataFrame in one pass; the old row-by-row `.loc` growth was O(n^2) over ~10^5 cells.
    all_concat = gpd.GeoDataFrame(geometry=gpd.GeoSeries(geoms, index=index))
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
        "Loaded elements: " + ", ".join(_element_names(sdata)),
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

# Rename transcript columns.
# In newer pandas/dask, dask's DataFrame.rename() no longer propagates the frame's `.attrs`, so it
# returns a Points frame stripped of `transform` + `spatialdata_attrs`. Reassigning that bare frame
# to sdata["transcripts"] then fails PointsModel validation with
#   "dask.dataframe.core.DataFrame.attrs does not contain `transform`" (scverse/spatialdata#996).
# We snapshot the valid `.attrs` before renaming and restore them afterwards, which reproduces the
# old (attrs-propagating) behaviour exactly and keeps the element model unchanged.
transcripts = sdata["transcripts"]
saved_attrs = deepcopy(dict(transcripts.attrs))

rename_map = {"global_z": "z", "transcript_id": "ensembl_id"}
renamed_gene = "gene" in transcripts.columns
if renamed_gene:
    # The feature_key column ('gene') is renamed too; update the feature_key attr to match,
    # otherwise validation KeyErrors on the now-missing 'gene' column.
    rename_map["gene"] = "feature_name"
    saved_attrs["spatialdata_attrs"]["feature_key"] = "feature_name"

renamed = transcripts.rename(columns=rename_map)
# Restore by mutating the frame's live `.attrs` dict in place (rather than reassigning `.attrs`),
# so this depends only on the getter returning the live dict -- the same assumption the rest of the
# loader already makes.
renamed.attrs.clear()
renamed.attrs.update(saved_attrs)
sdata["transcripts"] = renamed
print(
    datetime.now() - t0,
    "Renamed transcripts columns 'global_z' -> 'z', 'transcript_id' -> 'ensembl_id'"
    + (" and 'gene' -> 'feature_name'" if renamed_gene else ""),
    flush=True,
)

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

if n_iter == 1:
    # Common case (every current dataset has <65535 cells): a single rasterize call.
    # Keep the lazy DataArray exactly as before.
    print(datetime.now() - t0, "Rasterizing all cells in a single pass", flush=True)
    labels_image = sd.rasterize(
        sdata["cell_boundaries"], ["x", "y"], **rasterize_args
    )
else:
    # >65535 cells: rasterize in chunks and merge. Two subtleties make the naive in-place
    # merge wrong (it was silently broken before; this path had never run on real data since
    # no dataset exceeds 65535 cells):
    #   1. sd.rasterize returns a *lazy* (dask-backed) DataArray whose `.values` is a computed
    #      copy, so the old `labels_image.values[mask] = ...` write never persisted. We
    #      materialise each chunk and merge in a plain numpy array instead.
    #   2. return_regions_as_labels labels each chunk 1..len(chunk) *positionally*, so chunks
    #      would collide on the same label. We offset each chunk by its global start index
    #      (and use uint32, since labels exceed uint16 past 65535 cells).
    combined = None
    template = None
    for i in range(n_iter):
        start = i * N
        end = min((i + 1) * N, n_cells)
        print(
            datetime.now() - t0,
            f"Rasterizing iteration {i + 1}/{n_iter} (cells {start}..{end})",
            flush=True,
        )
        chunk = sd.rasterize(
            sdata["cell_boundaries"].iloc[start:end], ["x", "y"], **rasterize_args
        )
        chunk_np = np.asarray(chunk.data)
        if combined is None:
            combined = chunk_np.astype("uint32")
            template = chunk
        else:
            mask = chunk_np > 0
            combined[mask] = chunk_np[mask].astype("uint32") + start
    # Rebuild a dask-backed labels element (spatialdata rejects a numpy-backed one) while
    # preserving the template's transform and attrs (incl. label_index_to_category).
    labels_image = template.copy(
        data=da.from_array(combined, chunks=template.data.chunksize)
    )

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
