import spatialdata as sd
import spatialdata_plot as pl
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import cv2
import dask.dataframe as dd  #dask.dataframe should only be imported after spatial data to avoid prevent dependency conflicts
import scanpy as sc
import pandas as pd
import natsort
import os 
from bidcell import BIDCellModel


# preprocessing test data to get required input files for BIDcell
sdata = sd.read_zarr("dataset.zarr")
sdata_genes = sdata['transcripts']["feature_name"].unique().compute().sort_values().tolist()
# Extracting DAPI image from dataset.zarr
image_pyramid = []
img = sdata["morphology_mip"]["/scale0"]["image"].values  # Convert dask array to numpy
img = np.squeeze(img)  # Remove singleton channel dimension (c:1)
image_pyramid.append(img)

with tifffile.TiffWriter( "morphology_mip_pyramidal.tiff", bigtiff=True) as tiff:
    for img in image_pyramid:
        tiff.write(img, photometric="minisblack", resolution=(1, 1))


#Converting h5ad single cell reference to .csv
adata = sc.read_h5ad("dataset.h5ad")
shared_genes = [g for g in sdata_genes if g in adata.var["feature_name"].values] #crucial to avoid key error due to discrepencies in the genes present
adata = adata[:, adata.var["feature_name"].isin(shared_genes)]

# Set var_names to feature_name
adata.var_names = adata.var["feature_name"].astype(str)
sc_ref = pd.DataFrame(
    data=adata[:, shared_genes].layers["normalized"].toarray(), 
    columns=shared_genes, 
    index=range(adata.n_obs)
)
celltypes = adata.obs['cell_type'].unique().tolist()
cell_type_col = adata.obs['cell_type'].astype('category')
sc_ref["ct_idx"] = cell_type_col.cat.codes.values
sc_ref["cell_type"] = cell_type_col.values
sc_ref["atlas"] = "custom"
sc_ref.to_csv( "scref.csv")

# generating transcript map .csv from test data
transcript = sdata["transcripts"].compute()
transcript=pd.DataFrame(transcript)
#transcript.to_csv(filename="transcript.csv.gz", single_file=True, compression='gzip') #to generate transcript file without the filter
transcript[transcript["feature_name"].isin(shared_genes)].to_csv("transcript.csv.gz",compression='gzip'
)


#generating positive and negative marker files from the single cell reference
# defining the function generate_markers

def generate_markers(ref_df, max_overlaps_pos=4, max_overlaps_neg=15):
    """
    Generate positive and negative marker dataframes from reference data.
    
    Args:
        ref_df (pd.DataFrame): Reference dataframe with gene expression data and cell type info
        max_overlaps_pos (int): Maximum number of cell types that can share a positive marker
        max_overlaps_neg (int): Maximum number of cell types that can share a negative marker
    
    Returns:
        tuple: (df_pos, df_neg) - DataFrames containing positive and negative markers
    """

    n_genes = ref_df.shape[1] - 3
    cell_types = natsort.natsorted(list(set(ref_df["cell_type"].tolist())))
    n_cell_types = len(cell_types)
    
    ref_expr = ref_df.iloc[:, :n_genes].to_numpy()
    gene_names = ref_df.columns[:n_genes]
    ct_idx = ref_df["ct_idx"].to_numpy()

    # Generate negative markers
    pct_10 = np.percentile(ref_expr, 10, axis=1, keepdims=True)
    pct_10 = np.tile(pct_10, (1, n_genes))
    low_expr_true = np.zeros(pct_10.shape)
    low_expr_true[ref_expr <= pct_10] = 1

    low_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        low_expr_true_ct = low_expr_true[rows]
        low_expr_true_agg[ct, :] = np.prod(low_expr_true_ct, axis=0)

    overlaps = np.sum(low_expr_true_agg, 0)
    too_many = np.where(overlaps > max_overlaps_neg)[0]
    low_expr_true_agg[:, too_many] = 0
    df_neg = pd.DataFrame(low_expr_true_agg, index=cell_types, columns=gene_names)

    # Generate positive markers
    pct_90 = np.percentile(ref_expr, 90, axis=1, keepdims=True)
    pct_90 = np.tile(pct_90, (1, n_genes))
    high_expr_true = np.zeros(pct_90.shape)
    high_expr_true[ref_expr >= pct_90] = 1

    high_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        high_expr_true_ct = high_expr_true[rows]
        high_expr_true_agg[ct, :] = np.prod(high_expr_true_ct, axis=0)

    overlaps = np.sum(high_expr_true_agg, 0)
    too_many = np.where(overlaps > max_overlaps_pos)[0]
    high_expr_true_agg[:, too_many] = 0
    df_pos = pd.DataFrame(high_expr_true_agg, index=cell_types, columns=gene_names)

    return df_pos, df_neg

# generate positive and negative marker files 
df_pos, df_neg = generate_markers(sc_ref, max_overlaps_pos=4, max_overlaps_neg=15)
df_pos.to_csv( "pos_marker.csv")
df_neg.to_csv("neg_marker.csv")


# Setting up and running BIDCell
model = BIDCellModel("testdata.yaml")
model.run_pipeline() 

