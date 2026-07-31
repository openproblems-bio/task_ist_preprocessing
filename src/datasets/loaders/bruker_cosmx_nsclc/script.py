

"""
Notes about supported data structure (CosMx NSCLC / lung cancer, "Version 3"):

The general `bruker_cosmx` loader handles the mouse-brain / liver CosMx releases (one big zip,
optionally a separate flat-files zip). The NSCLC lung-cancer release is shaped differently and
is handled here, but this loader now follows the SAME pattern as the general one: the raw
archives are mirrored to S3 once (see scripts/create_resources/spatial/mirror_bruker_to_s3.sh)
and streamed in place at load time — only the members sopa actually needs are extracted, so the
enormous per-z-plane morphology stacks never touch scratch. The archives are passed as strings
(not staged Viash files) so Nextflow does not copy them first.

Each NSCLC sample ships TWO gzipped tarballs on the NanoString share:

1. <sample>+SMI+Flat+data.tar.gz  (par["input_raw"]) --> decompressed:

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

2. <sample>+RawMorphologyImages.tar.gz  (par["input_morphology"]) --> decompressed:

    └── <sample> 2/<sample>-RawMorphologyImages ( = morphology dir)
        ├── <some_id>_F001_Z001.TIF
        ├── ...
        ├── <some_id>_F001_Z008.TIF
        ├── ...
        └── <some_id>_F<last_fov_id>_Z008.TIF

The morphology images have one TIF per z-plane; sopa expects a single 2D image per FOV under
`DATA_DIR / "Morphology2D"`. We keep only the 4th plane (Z004) — this is the important filter,
it drops ~7/8 of the (very large) morphology archive — and rename e.g.
`<some_id>_F001_Z004.TIF` to `<some_id>_F001.tif` under `Morphology2D`.

The NSCLC flat files carry no `-polygons.csv`, so `cell_boundaries` is derived from the stitched
label image (`sd.to_polygons`), matching the original NSCLC loader behaviour.
"""

import os
import re
import shutil
import tarfile
from pathlib import Path
from datetime import datetime
import spatialdata as sd
import sopa

## VIASH START
par = {
    # Mirrored from
    # https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung9_Rep1/Lung9_Rep1+SMI+Flat+data.tar.gz
    "input_raw": "s3://openproblems-data/resources/raw_data/bruker_cosmx/Lung9_Rep1+SMI+Flat+data.tar.gz",
    # https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung9_Rep1/Lung9_Rep1+RawMorphologyImages.tar.gz
    "input_morphology": "s3://openproblems-data/resources/raw_data/bruker_cosmx/Lung9_Rep1+RawMorphologyImages.tar.gz",
    "segmentation_id": ["cell"],
    "output": "output.zarr",
    "dataset_id": "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung9_rep1",
    "dataset_name": "value",
    "dataset_url": "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "human",
}
meta = {
    "temp_dir": "./temp/datasets/bruker_cosmx_nsclc/temp",
}
## VIASH END

assert ("cell" in par["segmentation_id"]) and (len(par["segmentation_id"]) == 1), "Currently cell labels are definitely assigned in this script. And cosmx does not provide other segmentations."

t0 = datetime.now()

def log(msg):
    print(datetime.now() - t0, msg, flush=True)

# Define temp dir
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)

# Only the 4th morphology z-plane is kept (see module docstring). The regex is
# case-insensitive because the raw files use `.TIF`.
_Z_PLANE_RE = re.compile(r"_Z004\.tif$", re.I)


def _open_source(src):
    """Seekable binary handle for a local path or an s3:// URL (streamed, never fully downloaded)."""
    s = str(src)
    if s.startswith("s3://"):
        import s3fs
        # openproblems-data is public; anonymous read is deterministic and needs no creds.
        fs = s3fs.S3FileSystem(anon=True, default_block_size=32 * 2**20)
        return fs.open(s, "rb")
    return open(s, "rb")


def _member_wanted_default(name: str) -> bool:
    parts = name.split("/")
    return not (name.startswith("__MACOSX/") or parts[-1] == ".DS_Store")


def _member_wanted_morphology(name: str) -> bool:
    return _member_wanted_default(name) and bool(_Z_PLANE_RE.search(name))


def extract_targz(input_targz, output_dir: Path, wanted=_member_wanted_default):
    """
    Stream a .tar.gz (local path or s3:// URL) into output_dir, writing only the members for
    which `wanted(name)` is True. The archive is read once, sequentially (gzip is not seekable),
    but only wanted members are decoded to disk — so for the morphology tarball only the ~1/8 of
    files that are Z004 planes ever hit scratch. The archive itself is never fully downloaded to
    a local file first.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    kept = skipped = 0
    kept_bytes = 0
    with _open_source(input_targz) as fh, tarfile.open(fileobj=fh, mode="r|gz") as tar:
        for member in tar:
            if not member.isfile() or not wanted(member.name):
                skipped += 1
                continue
            # Flatten: strip any leading directory components, keep a stable relative path so
            # globbing below is layout-independent.
            target = output_dir.joinpath(*Path(member.name).parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            if extracted is None:
                skipped += 1
                continue
            with extracted as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
            kept += 1
            kept_bytes += member.size
    log(f"Extracted {kept} files ({kept_bytes / 1e9:.2f} GB); skipped {skipped} members")


#########################################
# Extract raw + morphology archives     #
#########################################

log("Extract tar.gz of flat files and cell labels")
FLAT_EXTRACTED = TMP_DIR / "input_raw"
extract_targz(par["input_raw"], FLAT_EXTRACTED)

# DATA_DIR is the directory that holds the flat CSVs + CellLabels (the "*-Flat_files_and_images"
# folder). Locate it via the tx flat file rather than hard-coding the sample-name nesting.
tx_files = sorted(FLAT_EXTRACTED.rglob("*_tx_file.csv"))
assert tx_files, (
    f"No '*_tx_file.csv' found under {FLAT_EXTRACTED}. "
    f"Contents: {[p.name for p in FLAT_EXTRACTED.rglob('*')][:50]}"
)
DATA_DIR = tx_files[0].parent
log(f"Detected flat-files directory: {DATA_DIR}")

# Check that all flat files are present
FLAT_FILES_ENDINGS = ["_exprMat_file.csv", "_fov_positions_file.csv", "_metadata_file.csv", "_tx_file.csv"]
for ending in FLAT_FILES_ENDINGS:
    present = any(f.name.endswith(ending) for f in DATA_DIR.iterdir())
    log(f"Flat file with ending {ending} is {'present' if present else 'MISSING'}")

log("Extract tar.gz of morphology images (Z004 plane only)")
MORPH_EXTRACTED = TMP_DIR / "input_morphology"
extract_targz(par["input_morphology"], MORPH_EXTRACTED, wanted=_member_wanted_morphology)

# Move the Z004 morphology tifs into DATA_DIR / "Morphology2D", renamed so sopa reads one 2D
# image per FOV (e.g. <id>_F001_Z004.TIF -> <id>_F001.tif).
log(f"Move Z004 morphology images to {DATA_DIR / 'Morphology2D'} and rename")
(DATA_DIR / "Morphology2D").mkdir(parents=True, exist_ok=True)
z004_tifs = sorted(MORPH_EXTRACTED.rglob("*_Z004.[Tt][Ii][Ff]"))
assert z004_tifs, f"No '*_Z004.TIF' morphology images found under {MORPH_EXTRACTED}"
for image_file in z004_tifs:
    dest_name = _Z_PLANE_RE.sub(".tif", image_file.name)
    shutil.move(str(image_file), DATA_DIR / "Morphology2D" / dest_name)
log(f"Moved {len(z004_tifs)} morphology images")

# sopa expects a CellLabels/ folder. The NSCLC flat archive ships one directly, but some CosMx
# exports only ship per-FOV label tifs inside FOV* folders — handle both, matching the general
# bruker_cosmx loader. See https://github.com/gustaveroussy/sopa/issues/285
labels_dir = DATA_DIR / "CellLabels"
fov_label_tifs = sorted(DATA_DIR.glob("FOV*/CellLabels_F*.tif"))
if labels_dir.is_dir():
    log("CellLabels folder already present in raw data")
elif fov_label_tifs:
    log(f"Symlink {len(fov_label_tifs)} per-FOV CellLabels tifs into {labels_dir}")
    labels_dir.mkdir(parents=True, exist_ok=True)
    for tif in fov_label_tifs:
        os.symlink(tif.resolve(), labels_dir / tif.name)
else:
    raise AssertionError(f"No CellLabels folder and no per-FOV CellLabels tifs found in {DATA_DIR}")

#############################################
# Memory optimization: lean transcript read #
#############################################
# sopa's CosMx reader loads the ENTIRE `*_tx_file.csv` into a pandas DataFrame in one
# `pd.read_csv` (sopa/io/reader/cosmx.py::_CosMXReader.read_transcripts). The string columns
# ('target', 'CellComp', ...) come in as object dtype (~60 B/row) — the peak of the loader.
# We wrap the `read_csv` used inside sopa's cosmx module so that, *only for the transcripts
# file*, we (a) read the 'target' feature column as categorical and (b) keep just the columns
# sopa's stitched-FOV path (read_transcripts + _get_global_cell_id) and the downstream schema
# actually use. Every other sopa read is left untouched. (Same shim as the general loader.)
from sopa.io.reader import cosmx as _cosmx_mod

_real_pd = _cosmx_mod.pd
_real_read_csv = _real_pd.read_csv

_TX_KEEP = ["fov", "cell_ID", "target", "x_global_px", "y_global_px", "z"]
_TX_REQUIRED = {"target", "x_global_px", "y_global_px", "fov", "cell_ID"}


def _lean_read_csv(filepath_or_buffer, *args, **kwargs):
    name = str(filepath_or_buffer)
    if "_tx_file.csv" in name and "usecols" not in kwargs:
        header = _real_read_csv(
            filepath_or_buffer, nrows=0,
            compression=kwargs.get("compression", "infer"),
        )
        cols = set(header.columns)
        if _TX_REQUIRED.issubset(cols):
            kwargs["usecols"] = [c for c in _TX_KEEP if c in cols]
            dtype = dict(kwargs.get("dtype") or {})
            dtype.setdefault("target", "category")
            kwargs["dtype"] = dtype
            log(f"Lean transcript read of {Path(name).name}: usecols={kwargs['usecols']}, target=category")
    return _real_read_csv(filepath_or_buffer, *args, **kwargs)


class _PandasReadCsvProxy:
    """Forwards every attribute to pandas, but leans the transcripts read_csv."""
    def __init__(self, real):
        self.__dict__["_real"] = real
    def read_csv(self, *args, **kwargs):
        return _lean_read_csv(*args, **kwargs)
    def __getattr__(self, item):
        return getattr(self._real, item)


_cosmx_mod.pd = _PandasReadCsvProxy(_real_pd)

#########################################
# Convert raw files to spatialdata zarr #
#########################################

log("Convert raw files to spatialdata zarr")

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
    cells_polygons=False,  # NSCLC flat files ship no -polygons.csv; derived from labels below
    flip_image=False,
)

# Retrieve polygons from segmentation (no polygons file in flat files present for NSCLC samples)
log("Derive cell polygons from stitched labels")
sdata["cells_polygons"] = sd.to_polygons(sdata["stitched_labels"])

###############
# Rename keys #
###############
log("Rename keys")
elements_renaming_map = {
    "stitched_image"     : "morphology_mip",
    "stitched_labels"    : "cell_labels",
    "points"             : "transcripts",
    "cells_polygons"     : "cell_boundaries",
    #"table"              : "metadata",
}

for old_key, new_key in elements_renaming_map.items():
    if old_key not in sdata:
        log(f"Element '{old_key}' not present, skipping rename to '{new_key}'")
        continue
    sdata[new_key] = sdata[old_key]
    del sdata[old_key]

# Add transcript columns with the names expected downstream.
sdata['transcripts']["cell_id"] = sdata['transcripts']["global_cell_id"]
sdata['transcripts']["feature_name"] = sdata['transcripts']["target"]

#########################################
# Throw out all channels except of DAPI #
#########################################
log("Throw out all channels except of 'DAPI'")

# TODO: in the future we want to keep PolyT and Cellbound1/2/3 stains. Note however, that somehow saving or plotting the sdata fails when
#       these channels aren't excluded, not sure why...
sdata["morphology_mip"] = sdata["morphology_mip"].sel(c=["DAPI"])

##############################
# Add info to metadata table #
##############################
log("Add metadata to table")
# The common iST API requires a per-cell unique `cell_id` in obs. sopa keys the table by the
# global cell id (the index), so expose that as `cell_id`. But the CosMx metadata already
# carries a per-FOV-*local* `cell_ID` column, and obs cannot hold two keys that differ only in
# case ('cell_ID' vs 'cell_id') — spatialdata's zarr writer rejects that as an invalid name.
# Rename the vendor local id out of the way first.
if "cell_ID" in sdata["table"].obs.columns:
    sdata["table"].obs.rename(columns={"cell_ID": "cell_ID_local"}, inplace=True)
sdata["table"].obs["cell_id"] = sdata["table"].obs.index.astype(str)
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism", "segmentation_id"]:
    sdata["table"].uns[key] = par[key]

#########
# Write #
#########
log(f"Writing to {par['output']}")

sdata.write(par["output"])

log(f"Done in {datetime.now() - t0}")
