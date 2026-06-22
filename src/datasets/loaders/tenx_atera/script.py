## code author: Florian Heyl
import shutil
import os
import zipfile
import tempfile
from pathlib import Path
from spatialdata_io import xenium

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
    "segmentation_id": par["segmentation_id"],
}
for key, value in new_uns.items():
    sdata.tables["table"].uns[key] = value

print(f"Output: {sdata}", flush=True)

print(f"Writing to '{par['output']}'", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])

sdata.write(par["output"])

print("Done", flush=True)