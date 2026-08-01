import spatialdata as sd
import sys
import numpy as np
import pandas as pd
import xarray as xr
import dask
import tacco
import anndata as ad
import argparse
from scipy.sparse import csr_matrix


def generate_adata(input_spots, cell_id_col="cell", gene_col="Gene"):
    """Aggregate per-transcript assignments into a cell x gene AnnData.

    Reimplementation of txsim.preprocessing.generate_adata: txsim's version is
    incompatible with anndata>=0.12 (it passes the removed `dtype=` argument to
    AnnData() and uses per-row item assignment `adata[cell, :] = ...`, which
    anndata no longer supports). This fills the count matrix directly and
    produces an identical output structure (X, layers['raw_counts'],
    obs/var stats, uns['spots'/'pct_noise']). Kept in sync with the copy in
    src/methods_count_aggregation/basic_count_aggregation/script.py.
    """
    spots = input_spots.copy()
    pct_noise = sum(spots[cell_id_col] <= 0) / len(spots[cell_id_col])
    spots_raw = spots.copy()  # kept in uns['spots']; 0 (background) -> None
    spots_raw.loc[spots_raw[cell_id_col] == 0, cell_id_col] = None
    spots = spots[spots[cell_id_col] > 0]

    cell_ids = pd.unique(spots[cell_id_col])
    genes = pd.unique(spots[gene_col])
    cell_pos = {c: i for i, c in enumerate(cell_ids)}

    # Populate the count matrix + centroids per cell (no AnnData item assignment).
    # feature_name may be categorical, so value_counts() can return unobserved
    # categories; reindex to `genes` (fill 0) to align to the var order, mirroring
    # txsim's `.reindex(var_names, fill_value=0)`.
    X = np.zeros((len(cell_ids), len(genes)), dtype=np.float32)
    centroid_x = np.zeros(len(cell_ids))
    centroid_y = np.zeros(len(cell_ids))
    for cell_id, grp in spots.groupby(cell_id_col, sort=False):
        row = cell_pos[cell_id]
        cts = grp[gene_col].value_counts().reindex(genes, fill_value=0)
        X[row, :] = cts.values
        centroid_x[row] = grp["x"].mean()
        centroid_y[row] = grp["y"].mean()

    adata = ad.AnnData(csr_matrix(X))
    adata.obs["cell_id"] = cell_ids
    adata.obs_names = [f"{i:d}" for i in cell_ids]
    adata.var_names = genes
    adata.obs["centroid_x"] = centroid_x
    adata.obs["centroid_y"] = centroid_y

    adata.uns["spots"] = spots_raw
    adata.uns["pct_noise"] = pct_noise
    adata.layers["raw_counts"] = adata.X.copy()

    adata.obs["n_counts"] = np.ravel(adata.layers["raw_counts"].sum(axis=1))
    adata.obs["n_genes"] = adata.layers["raw_counts"].getnnz(axis=1)
    adata.var["n_counts"] = np.ravel(adata.layers["raw_counts"].sum(axis=0))
    adata.var["n_cells"] = adata.layers["raw_counts"].getnnz(axis=0)
    return adata


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process input files and generate output files for FastReseg input"
    )
    parser.add_argument('input_path_ist', help='Path to the input file')
    parser.add_argument('input_segmentation_path', help='Path to the input segmentation file')
    parser.add_argument('input_sc_reference_path', help='Path to the input single-cell reference file')
    parser.add_argument('output_path_counts', help='Path for the output TSV file')
    parser.add_argument('output_path_transcripts', help='Path for the output TIF file')
    parser.add_argument('output_path_cell_type', help='Output cell type specification')
    
    return parser.parse_args()

### parsing arguments
args = parse_arguments()
print("args:")
print(args)
input_path = args.input_path_ist
input_segmentation_path = args.input_segmentation_path
input_sc_reference_path = args.input_sc_reference_path
print("path")
print(input_sc_reference_path)
output_path_counts = args.output_path_counts
output_path_transcripts = args.output_path_transcripts
output_path_cell_type = args.output_path_cell_type

## potential other parameters (TODO - make configurable)
sc_celltype_key = 'cell_type'

### reading the data in
sdata = sd.read_zarr(input_path)

# The transcripts points can carry a non-unique index: multi-partition parquet where
# each partition's RangeIndex restarts at 0, or concatenated "combined" replicates that
# each kept their own row index. dask-expr aligns on the index when the cell_id column is
# assigned below (`sdata['transcripts']["cell_id"] = ...`), and pandas cannot reindex a
# duplicate-labelled axis -> the following sd.transform then fails with
# "cannot reindex on an axis with duplicate labels". Reset to a unique RangeIndex eagerly
# here: a lazy .reset_index() does not help (the duplicate index already poisons the dask
# graph), and re-parsing via PointsModel.parse breaks the from_dask_array(index=...) assign
# below with a "Missing dependency" graph error, so materialise to pandas and re-wrap as a
# single from_pandas partition, re-attaching the coordinate transformations. No-op when the
# index is already unique.
if not sdata['transcripts'].index.compute().is_unique:
    print('Transcripts index has duplicate labels; resetting to a unique index', flush=True)
    _transforms = sd.transformations.get_transformation(sdata['transcripts'], get_all=True)
    _fixed = dask.dataframe.from_pandas(
        sdata['transcripts'].compute().reset_index(drop=True), npartitions=1
    )
    sd.transformations.set_transformation(_fixed, _transforms, set_all=True)
    sdata['transcripts'] = _fixed

### reading in basic segmentation
sdata_segm = sd.read_zarr(input_segmentation_path)
segmentation_coord_systems = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).keys()

# In case of a translation transformation of the segmentation (e.g. crop of the data), we need to adjust the transcript coordinates
trans = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True)['global'].inverse()

transcripts = sd.transform(sdata['transcripts'], to_coordinate_system='global')
transcripts = sd.transform(transcripts, trans, 'global')

print('Assigning transcripts to cell ids', flush=True)
y_coords = transcripts.y.compute().to_numpy(dtype=np.int64)
x_coords = transcripts.x.compute().to_numpy(dtype=np.int64)
if isinstance(sdata_segm["segmentation"], xr.DataTree):
    label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy()
else:
    label_image = sdata_segm["segmentation"].to_numpy()
# Guard against an empty segmentation. If the label image has no cells (all zeros)
# every transcript is assigned to background (cell_id 0), generate_adata drops them
# all, and the downstream tacco annotation then divides by a zero count norm and
# fails with a cryptic "divide by zero" instead of a useful message. Observed for the
# vizgen lung-cancer merscope dataset, whose polygon->label rasterization produced an
# all-zero image (for custom_segmentation this means the dataset's `cell_labels` is
# empty; re-process the dataset). Mirrors the guard in methods_transcript_assignment/pciseq.
if label_image.max() == 0:
    raise ValueError(
        f"Segmentation '{input_segmentation_path}' contains no cells: the label image "
        f"(shape={label_image.shape}, dtype={label_image.dtype}) is entirely zero. FastReseg "
        f"cannot assign transcripts to cells "
        f"(for custom_segmentation this means the dataset's `cell_labels` is empty)."
    )
# Clamp coords to the label-image bounds: transcripts at the crop boundary can
# round a few pixels past the raster edge, or the transcript field of view can
# extend beyond the segmented raster (real data), so label_image[y, x] would
# raise IndexError. Matches basic_transcript_assignment; edge/background at the border.
n_oob = int(np.count_nonzero((y_coords < 0) | (y_coords >= label_image.shape[0]) | (x_coords < 0) | (x_coords >= label_image.shape[1])))
y_coords = np.clip(y_coords, 0, label_image.shape[0] - 1)
x_coords = np.clip(x_coords, 0, label_image.shape[1] - 1)
print(f"Clamped {n_oob}/{len(x_coords)} transcripts outside the {label_image.shape[0]}x{label_image.shape[1]} label image to its edge", flush=True)
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        label_image[y_coords, x_coords], chunks=tuple(sdata['transcripts'].map_partitions(len).compute())
    ), 
    index=sdata['transcripts'].index
)
sdata['transcripts']["cell_id"] = cell_id_dask_series

### extracting transcript ids
print('Transforming transcripts coordinates', flush=True)
transcripts = sd.transform(sdata['transcripts'], to_coordinate_system='global')

# .copy() defensively: for a single-partition points object, transcripts.compute()
# can return a reference to the dask graph's backing frame, so the in-place rename
# below would otherwise corrupt it and any later compute of `transcripts` would yield
# renamed columns that no longer match the dask _meta -> "Metadata mismatch".
transcripts_df = transcripts.compute().copy()

# Some datasets carry no per-molecule transcript id: e.g. the Vizgen MERSCOPE loader
# renames the vendor's `transcript_id` (a gene/ensembl id) to `ensembl_id`, leaving the
# transcripts without a `transcript_id` column. FastReseg itself does not use one
# (transID_coln = NULL in script.R), but output.py needs a `transcript_id` in the final
# output, so synthesise a unique one from the row position when it is absent.
if 'transcript_id' not in transcripts_df.columns:
    print('No transcript_id column found; synthesising a unique transcript_id', flush=True)
    transcripts_df['transcript_id'] = np.arange(len(transcripts_df), dtype=np.int64)

transcripts_df.rename(columns = {'feature_name': 'target',
'transcript_id': 'UMI_transID', 'cell_id': 'UMI_cellID'}, inplace = True)

transcripts_df = transcripts_df.loc[:, ['target', 'x', 'y', 'z', 'UMI_transID', 'UMI_cellID']]
transcripts_df.to_csv(output_path_transcripts)


#### aggregating counts per transcript, based on 
df = sdata['transcripts'].compute()
df.feature_name = df.feature_name.astype(str)

adata_sp = generate_adata(df, cell_id_col='cell_id', gene_col='feature_name') #TODO: x and y refers to a specific coordinate system. Decide which space we want to use here. (probably should be handled in the previous assignment step)
adata_sp.layers['counts'] = adata_sp.layers['raw_counts']
del adata_sp.layers['raw_counts']
adata_sp.var["gene_name"] = adata_sp.var_names
print(adata_sp.var_names[1:10])

# currently the function also saves the transcripts in the adata object, but this is not necessary here
del adata_sp.uns['spots']
del adata_sp.uns['pct_noise']


count_df = pd.DataFrame(adata_sp.X.toarray(), 
                       index=adata_sp.obs_names, 
                       columns=adata_sp.var_names)
count_df.to_csv(output_path_counts)

#### run cell annotation with tacco
adata_sc = ad.read_h5ad(input_sc_reference_path)
print(adata_sc.var_names[1:10])

# tacco annotates the cell x gene count matrix directly against the reference, so
# (unlike the previous ssam path) it needs neither the transcript-level spots nor a
# gene subsetting step -- tacco intersects genes internally. It expects raw counts on
# .X for both spatial and reference, mirroring methods_cell_type_annotation/tacco.
adata_sp.X = adata_sp.layers['counts']
adata_sc.X = adata_sc.layers['counts']

print('Annotating cell types', flush=True)
cell_type_assignment = tacco.tl.annotate(
    adata=adata_sp,
    reference=adata_sc,
    annotation_key=sc_celltype_key,
)

# tacco returns an n_obs x n_celltypes proportion matrix (aligned to adata_sp.obs);
# take the argmax per cell for the hard label the FastReseg R step consumes as `clust`.
cell_types = cell_type_assignment.columns
highest_score_idx = np.argmax(cell_type_assignment.to_numpy(), axis=1)
adata_sp.obs[sc_celltype_key] = cell_types[highest_score_idx]

# Preserve the exact CSV shape the R step reads (read.csv(row.names=1) -> a single
# column of per-cell labels indexed by cell_id).
cell_type_df = adata_sp.obs[sc_celltype_key].astype(str)
cell_type_df.to_csv(output_path_cell_type, header=True)