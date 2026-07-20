## code author: Florian Heyl
import shutil
import os
import zipfile
import tempfile
from pathlib import Path
from spatialdata_io import xenium

try:
    from xarray import DataTree
except ImportError:  # older xarray ships DataTree as a separate package
    from datatree import DataTree


def rechunk_uniform(sdata, chunk_size=1024):
    """Rechunk all image/label rasters to a uniform (regular) chunk grid.

    xenium() can return multiscale rasters whose dask chunks are irregular.
    When written, an irregular grid is encoded as a *rectilinear* chunk grid,
    which spatialdata's ome-zarr reader cannot read back (it accesses the
    ``.chunks`` attribute, which is undefined for rectilinear grids). Forcing a
    uniform chunk size along the spatial dims guarantees a regular grid.
    """
    spatial = ("z", "y", "x")
    for group in (sdata.images, sdata.labels):
        for key in list(group.keys()):
            elem = group[key]
            # multiscale rasters are DataTree (indexable by "scale0");
            # single-scale rasters are a plain DataArray.
            ref = elem["scale0"] if isinstance(elem, DataTree) else elem
            coords = list(ref.coords.keys())
            strategy = {c: chunk_size for c in coords if c in spatial}
            if "c" in coords:  # keep the channel axis as a single chunk
                strategy["c"] = -1
            group[key] = elem.chunk(strategy)

## VIASH START
par = {
    "input": "temp/datasets/10x_atera/WTA_Preview_FFPE_Breast_Cancer_xe_outs.zip",
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
    "output": "temp/datasets/10x_atera/breast/breast.zarr"
}
meta = {
    "cpus": 1,
    "temp_dir": None,
}

## VIASH END


def extract_zip(input_zip: Path, output_dir: Path, strip_root: bool = False):
    output_dir = Path(output_dir)
    with zipfile.ZipFile(input_zip, 'r') as zip_ref:
        members = zip_ref.infolist()

        roots = {Path(m.filename).parts[0] for m in members if m.filename.strip("/") and not m.filename.startswith("__MACOSX/")}
        if not (strip_root and len(roots) == 1):
            zip_ref.extractall(output_dir)
            return

        for member in members:
            if member.filename.startswith("__MACOSX/"):
                continue
            parts = Path(member.filename).parts[1:]
            if not parts:
                continue
            target = output_dir.joinpath(*parts)
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with zip_ref.open(member) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)


TMP_DIR = Path(meta["temp_dir"] or tempfile.mkdtemp())
TMP_DIR.mkdir(parents=True, exist_ok=True)

print("Extract input zip", flush=True)
input_extracted = TMP_DIR / "input"
extract_zip(Path(par["input"]), input_extracted, strip_root=True)

print(f"Files in extracted dir: {os.listdir(input_extracted)}", flush=True)

# read the data
sdata = xenium(
    path=input_extracted,
    n_jobs=meta["cpus"] or 1,
    cells_boundaries=True,
    nucleus_boundaries=True,
    morphology_focus=True,
    cells_as_circles=False,
)

# Keep exactly one morphology raster named "morphology_mip"; process_dataset
# renames it to the API-required "image" element. Atera output mirrors Xenium
# Onboard Analysis v4, which ships only a multi-channel "morphology_focus"
# (channel 0 is DAPI, which the segmentation methods use via image[0]) and no
# separate "morphology_mip". Prefer the MIP if present, otherwise fall back to
# morphology_focus so the dataset always has an image.
if "morphology_mip" not in sdata.images and "morphology_focus" in sdata.images:
    sdata["morphology_mip"] = sdata["morphology_focus"]
if "morphology_focus" in sdata.images:
    del sdata.images["morphology_focus"]

# Log the morphology image channel names so it's visible in the run log which
# channel is which — segmentation uses channel 0 (expected to be DAPI).
try:
    _mip = sdata["morphology_mip"]
    try:
        _arr = _mip["scale0"].image  # multiscale raster
    except (TypeError, KeyError):
        _arr = _mip  # single-scale DataArray
    _channels = [str(c) for c in _arr.coords["c"].values] if "c" in _arr.coords else "<no channel axis>"
    print(f"morphology image channels (channel 0 used for segmentation): {_channels}", flush=True)
except Exception as e:
    print(f"Could not read morphology image channel names: {e}", flush=True)

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

# Force a regular chunk grid so the written store is readable by spatialdata's
# ome-zarr reader (see rechunk_uniform). Do NOT enable array.rectilinear_chunks:
# that only permits writing the unreadable rectilinear grids we are avoiding.
print("Rechunking rasters to a uniform chunk grid", flush=True)
rechunk_uniform(sdata)

print(f"Output: {sdata}", flush=True)

print(f"Writing to '{par['output']}'", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])

sdata.write(par["output"])

print("Done", flush=True)