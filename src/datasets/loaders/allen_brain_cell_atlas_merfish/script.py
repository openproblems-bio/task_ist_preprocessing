import math
import os
import re
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
import tifffile
from shapely.geometry import Polygon
import spatialdata as sd
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel

## VIASH START
par = {
    "mouse": "mouse1_coronal",
    "experiment_id": "220329_wb3_co1_1_3z19R_merfish1",
    "abca_version": "20231215",
    "segmentation_id": ["cell"],
    "output": "output.zarr",
    "dataset_id": "allen_brain_cell_atlas_merfish/mouse1_coronal/rep1",
    "dataset_name": "Allen Brain Cell Atlas MERFISH mouse1_coronal",
    "dataset_url": "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/",
    "dataset_reference": "@article{Yao2023, author={Yao, Zizhen and others}, title={A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain}, journal={Nature}, year={2023}}",
    "dataset_summary": "Brain-wide MERFISH spatial transcriptomics atlas of mouse brain.",
    "dataset_description": "Brain-wide MERFISH spatial transcriptomics data from the Zhuang lab. Four mice (2 coronal, 2 sagittal) imaged with ~1100 gene panel across the full brain volume.",
    "dataset_organism": "Mus musculus",
    "keep_files": True,
}
meta = {"temp_dir": "/tmp"}
## VIASH END

# ─── Constants ────────────────────────────────────────────────────────────────

MICRONS_PER_PIXEL = 0.109
FRAME_SIZE_PIXELS = 2048
PERCENT_SCALING = 0.3
OVERLAP_SIZE = 213

BASE_URL = "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
DATA_DOWNLOAD_LINKS = {
    "mouse1_coronal": BASE_URL + "coronal_1/",
    "mouse2_coronal": BASE_URL + "coronal_2/",
    "mouse3_sagittal": BASE_URL + "sagittal_1/",
    "mouse4_sagittal": BASE_URL + "sagittal_2/",
}
MOUSE_TO_ABCA_ID = {
    "mouse1_coronal": "Zhuang-ABCA-1",
    "mouse2_coronal": "Zhuang-ABCA-2",
    "mouse3_sagittal": "Zhuang-ABCA-3",
    "mouse4_sagittal": "Zhuang-ABCA-4",
}
# Known per-experiment image channel prefix corrections
DOWNLOAD_PREFIX_CORRECTIONS = {
    ("mouse1_coronal", "220422_wb3_co1_6B_6z18R_merfish6"): "Epi-750s6-635s6-545s1-473s6-408s6",
    ("mouse1_coronal", "220411_wb3_co1_5B_6z19R_merfish5"): "Epi-750s6-635s6-545s1-473s6-408s6",
    ("mouse1_coronal", "220329_wb3_co1_24_3z19R_merfish5_dryhyb"): "Epi-750s3-635s3-545s1-473s3-408s3",
    ("mouse2_coronal", "220506_wb3_co2_10_5z18R_merfish2"): "Epi-750s5-650s5-560s1-488s5-405s5",
}

t0 = datetime.now()
mouse = par["mouse"]
experiment_id = par["experiment_id"]

TMP_DIR = Path(meta.get("temp_dir") or "/tmp") / "allen_brain_merfish"
TMP_DIR.mkdir(parents=True, exist_ok=True)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def download_if_missing(url, path):
    if not Path(path).exists():
        print(datetime.now() - t0, f"Downloading {url}", flush=True)
        urllib.request.urlretrieve(url, path)


def check_url(url):
    try:
        return urllib.request.urlopen(url).getcode() == 200
    except (urllib.error.HTTPError, urllib.error.URLError):
        return False


def compute_weight(location, dimensions, percent_scaling):
    """Linear blending weight for FOV overlap. Adapted from Fiji Stitching plugin."""
    min_distance = 1
    for dim in range(len(location)):
        value = max(1, min(location[dim] + 1, (dimensions[dim] - 1) - location[dim] + 1))
        img_area_blend = round(percent_scaling * 0.5 * dimensions[dim])
        value = value / img_area_blend if value < img_area_blend else 1
        min_distance *= value
    if min_distance == 1:
        return 1
    if min_distance <= 0:
        return 0.0000001
    return (math.cos((1 - min_distance) * math.pi) + 1) / 2


# ─── Dax file reader (Zhuang lab custom format) ────────────────────────────────

class DaxReader:
    def __init__(self, filename):
        self.filename = filename
        inf_path = os.path.splitext(filename)[0] + ".inf"
        self.image_height = self.image_width = self.number_frames = None
        self.bigendian = 0

        with open(inf_path) as f:
            for line in f:
                m = re.match(r"frame dimensions = ([\d]+) x ([\d]+)", line)
                if m:
                    self.image_height, self.image_width = int(m.group(2)), int(m.group(1))
                m = re.match(r"number of frames = ([\d]+)", line)
                if m:
                    self.number_frames = int(m.group(1))
                m = re.search(r" (big|little) endian", line)
                if m:
                    self.bigendian = int(m.group(1) == "big")

        if not self.image_height:
            self.image_height = self.image_width = 256

        self.fileptr = open(filename, "rb")

    def load_frame(self, frame_number):
        self.fileptr.seek(frame_number * self.image_height * self.image_width * 2)
        img = np.fromfile(self.fileptr, dtype="uint16", count=self.image_height * self.image_width)
        img = img.reshape(self.image_height, self.image_width)
        if self.bigendian:
            img.byteswap(True)
        return img

    def close(self):
        self.fileptr.close()


def download_fov_image(fov_idx, n_format, prefix, suffix, base_url, frames, stain):
    """Download one FOV dax file and return a (n_z, px, py) uint16 array."""
    fmt = f"{fov_idx:03d}" if n_format == 3 else f"{fov_idx:04d}"
    stem = prefix + fmt + suffix
    dax_path = TMP_DIR / (stem + ".dax")
    inf_path = TMP_DIR / (stem + ".inf")
    tif_path = TMP_DIR / (stem + f"_{stain}.tif")

    if tif_path.exists():
        return tifffile.imread(tif_path)

    if not dax_path.exists():
        urllib.request.urlretrieve(base_url + stem + ".dax", dax_path)
        urllib.request.urlretrieve(base_url + stem + ".inf", inf_path)

    reader = DaxReader(str(dax_path))
    img = np.array([reader.load_frame(f) for f in frames])
    reader.close()
    tifffile.imsave(tif_path, img)
    return img


def delete_dax_files():
    for f in TMP_DIR.glob("*.dax"):
        f.unlink()
    for f in TMP_DIR.glob("*.inf"):
        f.unlink()


# ─── Step 1: Experiment metadata ──────────────────────────────────────────────

print(datetime.now() - t0, "Loading experiment metadata...", flush=True)

exp_metadata_path = TMP_DIR / "experiment_metadata.csv"
download_if_missing(BASE_URL + "additional_files/experiment_metadata.csv", exp_metadata_path)

exp_metadata = pd.read_csv(exp_metadata_path)
exp_metadata = exp_metadata[exp_metadata["experiment_id"] != "220514_sa_8_merfish4_adaptor"]
exp_metadata.loc[
    exp_metadata["experiment_id"] == "220627_wb3_sa1_5_5z18R_merfish5", "experiment_id"
] = "220627_wb3_sa1_6_5z18R_merfish6"

exp_row = exp_metadata[
    (exp_metadata["animal_id"] == mouse) & (exp_metadata["experiment_id"] == experiment_id)
].iloc[0]
sample_id = exp_row["sample_id"]

# Load data organization file to find DAPI channel and z-plane indices
data_org_path = TMP_DIR / f"dataorganization_{experiment_id}.csv"
download_if_missing(
    BASE_URL + f"additional_files/dataorganization/{mouse}/dataorganization_{experiment_id}.csv",
    data_org_path,
)
data_org = pd.read_csv(data_org_path)
dapi_row = data_org[data_org["channelName"] == "DAPI"].iloc[0]
dapi_frames = [int(n) for n in re.findall(r"\d+", dapi_row["frame"])]
dapi_prefix = dapi_row["imageType"]
n_z_planes = len(dapi_frames)

# Apply known per-experiment prefix corrections
correction = DOWNLOAD_PREFIX_CORRECTIONS.get((mouse, experiment_id))
if correction:
    dapi_prefix = correction

# Determine FOV file naming convention (3 or 4-digit index, optional _0 suffix)
image_base_url = DATA_DOWNLOAD_LINKS[mouse] + experiment_id + "/data/"
n_format, download_suffix = 3, ""
test_stem = dapi_prefix + "_"
if not check_url(image_base_url + test_stem + "002.dax"):
    if check_url(image_base_url + test_stem + "002_0.dax"):
        download_suffix = "_0"
    elif check_url(image_base_url + test_stem + "0002.dax"):
        n_format = 4
    elif check_url(image_base_url + test_stem + "0002_0.dax"):
        n_format, download_suffix = 4, "_0"

dapi_file_prefix = dapi_prefix + "_"

# ─── Step 2: FOV layout from raw count matrix ─────────────────────────────────

print(datetime.now() - t0, "Loading raw count matrix...", flush=True)

matrix_path = TMP_DIR / f"{mouse}_raw_matrix.h5ad"
download_if_missing(
    BASE_URL + f"processed_data/counts_updated/raw_counts_{mouse}.h5ad",
    matrix_path,
)

raw_adata = anndata.read_h5ad(matrix_path)
sample_adata = raw_adata[raw_adata.obs["sample_id"] == sample_id].copy()

fov_df = sample_adata.obs[["fov", "fov_x", "fov_y"]].drop_duplicates().copy()
n_fovs = int(sample_adata.obs["fov"].max())
present_fovs = set(fov_df["fov"].values)
missing_fovs = {i for i in range(n_fovs) if i not in present_fovs}
print(datetime.now() - t0, f"FOVs: {n_fovs} total, {len(missing_fovs)} missing", flush=True)

frame_size_um = FRAME_SIZE_PIXELS * MICRONS_PER_PIXEL
xmin, xmax = fov_df["fov_x"].min(), fov_df["fov_x"].max()
ymin, ymax = fov_df["fov_y"].min(), fov_df["fov_y"].max()

img_w = int(((xmax - xmin) + frame_size_um) / MICRONS_PER_PIXEL)
img_h = int(((ymax - ymin) + frame_size_um) / MICRONS_PER_PIXEL)


def um_to_px_x(x):
    # x-axis is flipped in the stitched image
    return img_w - int((x - xmin) / MICRONS_PER_PIXEL)


def um_to_px_y(y):
    return int((y - ymin) / MICRONS_PER_PIXEL)


fov_df["x_min"] = fov_df["fov_x"].apply(um_to_px_x)
fov_df["x_max"] = fov_df["x_min"] + FRAME_SIZE_PIXELS
fov_df["y_min"] = fov_df["fov_y"].apply(um_to_px_y)
fov_df["y_max"] = fov_df["y_min"] + FRAME_SIZE_PIXELS

# ─── Step 3: Stitch DAPI image ────────────────────────────────────────────────

print(datetime.now() - t0, "Stitching DAPI image...", flush=True)

# Precompute blending weight matrices for overlap regions
w1 = np.zeros((OVERLAP_SIZE, FRAME_SIZE_PIXELS))
w2 = np.zeros((OVERLAP_SIZE, FRAME_SIZE_PIXELS))
for xi in range(FRAME_SIZE_PIXELS):
    for yi in range(FRAME_SIZE_PIXELS):
        if xi >= FRAME_SIZE_PIXELS - OVERLAP_SIZE:
            w1[xi - (FRAME_SIZE_PIXELS - OVERLAP_SIZE), yi] = max(
                0.00001, compute_weight([xi, yi], [FRAME_SIZE_PIXELS, FRAME_SIZE_PIXELS], PERCENT_SCALING)
            )
        if xi < OVERLAP_SIZE:
            w2[xi, yi] = max(
                0.00001, compute_weight([xi, yi], [FRAME_SIZE_PIXELS, FRAME_SIZE_PIXELS], PERCENT_SCALING)
            )

complete_img = np.zeros((n_z_planes, img_w, img_h), dtype="uint16")

for fov in range(n_fovs):
    if fov in missing_fovs:
        continue

    row = fov_df[fov_df["fov"] == fov].iloc[0]
    xp, xp_max = int(row["x_min"]), int(row["x_max"])
    yp, yp_max = int(row["y_min"]), int(row["y_max"])

    fov_img = download_fov_image(fov, n_format, dapi_file_prefix, download_suffix, image_base_url, dapi_frames, "DAPI")

    # Find neighbors for linear blending at overlaps
    fov_right = fov_df[
        (fov_df["x_min"] == xp) & (fov_df["y_min"] > yp) & (fov_df["y_min"] < yp_max)
    ]["fov"].values
    fov_bottom = fov_df[
        (fov_df["y_min"] == yp) & (fov_df["x_min"] > xp) & (fov_df["x_min"] < xp_max)
    ]["fov"].values
    fov_br = fov_df[
        (fov_df["y_min"] > yp) & (fov_df["y_min"] < yp_max) &
        (fov_df["x_min"] > xp) & (fov_df["x_min"] < xp_max)
    ]["fov"].values

    right = len(fov_right) > 0 and fov_right[0] not in missing_fovs
    bottom = len(fov_bottom) > 0 and fov_bottom[0] not in missing_fovs

    if bottom:
        bot_img = download_fov_image(int(fov_bottom[0]), n_format, dapi_file_prefix, download_suffix, image_base_url, dapi_frames, "DAPI")
        fov_img[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:, :] = (
            np.multiply(fov_img[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:, :], w1) +
            np.multiply(bot_img[:, :OVERLAP_SIZE, :], w2)
        ) / (w1 + w2)

    if right:
        right_img = download_fov_image(int(fov_right[0]), n_format, dapi_file_prefix, download_suffix, image_base_url, dapi_frames, "DAPI")
        fov_img[:, :, FRAME_SIZE_PIXELS - OVERLAP_SIZE:] = (
            np.multiply(fov_img[:, :, FRAME_SIZE_PIXELS - OVERLAP_SIZE:], w1.T) +
            np.multiply(right_img[:, :, :OVERLAP_SIZE], w2.T)
        ) / (w1.T + w2.T)

    if len(fov_br) > 0 and fov_br[0] not in missing_fovs:
        br_img = download_fov_image(int(fov_br[0]), n_format, dapi_file_prefix, download_suffix, image_base_url, dapi_frames, "DAPI")
        sl_x = slice(FRAME_SIZE_PIXELS - OVERLAP_SIZE, None)
        sl_y = slice(FRAME_SIZE_PIXELS - OVERLAP_SIZE, None)
        if bottom and right:
            fov_img[:, sl_x, sl_y] = np.divide(
                np.multiply(fov_img[:, sl_x, sl_y], w1[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:]) +
                np.multiply(right_img[:, sl_x, :OVERLAP_SIZE], w1[:, :OVERLAP_SIZE]) +
                np.multiply(bot_img[:, :OVERLAP_SIZE, sl_y], w2[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:]) +
                np.multiply(br_img[:, :OVERLAP_SIZE, :OVERLAP_SIZE], w2[:, :OVERLAP_SIZE]),
                w1[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:] + w1[:, :OVERLAP_SIZE] +
                w2[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:] + w2[:, :OVERLAP_SIZE],
            )
        elif bottom:
            fov_img[:, sl_x, sl_y] = np.divide(
                np.multiply(fov_img[:, sl_x, sl_y], w1[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:]) +
                np.multiply(bot_img[:, :OVERLAP_SIZE, sl_y], w2[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:]) +
                np.multiply(br_img[:, :OVERLAP_SIZE, :OVERLAP_SIZE], w2[:, :OVERLAP_SIZE]),
                w1[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:] + w2[:, FRAME_SIZE_PIXELS - OVERLAP_SIZE:] + w2[:, :OVERLAP_SIZE],
            )

    # Skip overlap region already filled by a previously placed FOV
    if np.any(complete_img[:, xp:xp_max, yp:yp + OVERLAP_SIZE] != 0):
        fov_img = fov_img[:, :, OVERLAP_SIZE:]
        yp += OVERLAP_SIZE
    if np.any(complete_img[:, xp:xp + OVERLAP_SIZE, yp:yp_max - OVERLAP_SIZE] != 0):
        fov_img = fov_img[:, OVERLAP_SIZE:, :]
        xp += OVERLAP_SIZE

    complete_img[:, xp:xp + fov_img.shape[1], yp:yp + fov_img.shape[2]] = fov_img

    if fov % 10 == 0:
        delete_dax_files()
        print(datetime.now() - t0, f"Stitched {fov}/{n_fovs} FOVs", flush=True)

print(datetime.now() - t0, "Image stitching complete", flush=True)

dapi_mip = complete_img.max(axis=0)  # (img_w, img_h) max-projection over z
del complete_img

# ─── Step 4: Load and pixel-register transcripts ──────────────────────────────

print(datetime.now() - t0, "Loading transcripts...", flush=True)

spots_path = TMP_DIR / f"spots_{experiment_id}.csv"
download_if_missing(
    BASE_URL + f"processed_data/decoded_spots/{mouse}/spots_{experiment_id}.csv",
    spots_path,
)
spots_df = pd.read_csv(spots_path)
spots_df["x"] = spots_df["global_x"].apply(um_to_px_x).astype(float)
spots_df["y"] = spots_df["global_y"].apply(um_to_px_y).astype(float)
spots_df = spots_df.rename(columns={
    "target_gene": "feature_name",
    "barcode_id": "transcript_id",
    "global_z": "z",
})

# ─── Step 5: Load and pixel-register cell boundaries ─────────────────────────

print(datetime.now() - t0, "Loading cell boundaries...", flush=True)

boundaries_path = TMP_DIR / f"{experiment_id}_boundaries.csv"
download_if_missing(
    BASE_URL + f"processed_data/cell_boundaries_updated/{mouse}/{experiment_id}.csv",
    boundaries_path,
)
cell_df = pd.read_csv(boundaries_path)

def parse_polygon_z0(row):
    if not isinstance(row.get("boundaryX_z0"), str):
        return None
    xs = [float(v) for v in row["boundaryX_z0"].split(",")]
    ys = [float(v) for v in row["boundaryY_z0"].split(",")]
    return Polygon(zip([um_to_px_x(x) for x in xs], [um_to_px_y(y) for y in ys]))

print(datetime.now() - t0, "Parsing cell boundary polygons (z=0)...", flush=True)
cell_df["geometry"] = cell_df.apply(parse_polygon_z0, axis=1)
cell_df = cell_df.dropna(subset=["geometry"])
gdf = gpd.GeoDataFrame(cell_df[["Unnamed: 0", "geometry"]], geometry="geometry")
gdf = gdf.rename(columns={"Unnamed: 0": "cell_id"}).set_index("cell_id")
gdf.index = gdf.index.astype(str)

# ─── Step 6: Assign transcripts to cells (spatial join) ───────────────────────

print(datetime.now() - t0, "Assigning transcripts to cells...", flush=True)

spots_gdf = gpd.GeoDataFrame(
    spots_df,
    geometry=gpd.points_from_xy(spots_df["x"], spots_df["y"]),
)
joined = gpd.sjoin(
    spots_gdf,
    gdf.reset_index()[["cell_id", "geometry"]],
    how="left",
    predicate="within",
)
spots_df["cell_id"] = joined["cell_id"].values
print(
    datetime.now() - t0,
    f"Transcripts assigned: {spots_df['cell_id'].notna().sum()} / {len(spots_df)}",
    flush=True,
)

# ─── Step 7: Build metadata AnnData table ─────────────────────────────────────

print(datetime.now() - t0, "Loading cell type metadata...", flush=True)

abca_id = MOUSE_TO_ABCA_ID[mouse]
cell_meta_path = TMP_DIR / f"{mouse}_cell_metadata.csv"
download_if_missing(
    f"https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/{abca_id}/{par['abca_version']}/views/cell_metadata_with_cluster_annotation.csv",
    cell_meta_path,
)
cell_meta = pd.read_csv(cell_meta_path, index_col=0)
cell_meta.index = cell_meta.index.astype(str)

print(datetime.now() - t0, "Building metadata AnnData table...", flush=True)

# Keep only cells present in both the count matrix and the boundary polygons
gdf_cell_ids = set(gdf.index)
adata = sample_adata[sample_adata.obs.index.isin(gdf_cell_ids)].copy()
adata.obs.index = adata.obs.index.astype(str)

# Add cell type annotations
annotation_cols = ["cluster_alias", "class", "subclass", "supertype", "cluster"]
adata.obs = adata.obs.join(cell_meta[annotation_cols], how="left")

# Add spatial coordinates from cell boundary centroids
centroids = {str(k): (v.centroid.x, v.centroid.y) for k, v in gdf.geometry.items()}
adata.obsm["spatial"] = np.array([
    centroids.get(str(idx), (np.nan, np.nan)) for idx in adata.obs.index
])

adata.obs["cell_id"] = adata.obs.index.astype(str)
adata.obs["region"] = "cell_boundaries"

adata.var["gene_ids"] = adata.var.index.astype(str)
adata.var["feature_types"] = "Gene Expression"

for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference",
            "dataset_summary", "dataset_description", "dataset_organism", "segmentation_id"]:
    adata.uns[key] = par[key]

# ─── Step 8: Build SpatialData object ─────────────────────────────────────────

print(datetime.now() - t0, "Building SpatialData object...", flush=True)

identity = sd.transformations.Identity()

# Image: (c, x, y) — DAPI max-projection over z-planes
image = Image2DModel.parse(
    dapi_mip[np.newaxis, :, :],
    dims=("c", "x", "y"),
    c_coords=["DAPI"],
    transformations={"global": identity},
)

# Transcripts
transcript_cols = [c for c in ["x", "y", "z", "feature_name", "transcript_id", "cell_id"] if c in spots_df.columns]
transcripts = PointsModel.parse(
    spots_df[transcript_cols],
    coordinates={"x": "x", "y": "y", "z": "z"} if "z" in spots_df.columns else {"x": "x", "y": "y"},
    feature_key="feature_name",
    instance_key="transcript_id",
    transformations={"global": identity},
)

# Cell boundary polygons
shapes = ShapesModel.parse(
    gdf,
    transformations={"global": identity},
)

# Cell x gene table annotated to cell_boundaries
table = TableModel.parse(
    adata,
    region="cell_boundaries",
    region_key="region",
    instance_key="cell_id",
)

sdata = sd.SpatialData(
    images={"image": image},
    points={"transcripts": transcripts},
    shapes={"cell_boundaries": shapes},
    tables={"table": table},
)

# ─── Write ────────────────────────────────────────────────────────────────────

print(datetime.now() - t0, f"Writing to {par['output']}...", flush=True)

sdata.write(par["output"])

if not par.get("keep_files"):
    import shutil
    shutil.rmtree(TMP_DIR)

print(datetime.now() - t0, "Done", flush=True)
