import numpy as np
import anndata as ad
import spatialdata as sd
import os
import shutil

### VIASH START
par = {
  "input_sc": "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "input_sp": "resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
  "output_scrnaseq": "resources_test/task_ist_preprocessing/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad",
  "output_ist": "resources_test/task_ist_preprocessing/2023_10x_mouse_brain_xenium_rep1/dataset.zarr"
}
### VIASH END


def get_crop_coords(sdata, max_n_pixels=20000*20000): #50000*50000): 
    """Get the crop coordinates to subset the sdata to max_n_pixels
    
    Arguments
    ---------
    sdata: spatialdata.SpatialData
        The spatial data to crop
    max_n_pixels: int
        The maximum number of pixels to keep
        
    Returns
    -------
    crop_coords: dict
        The crop coordinates
    """
    
    _, h, w = sdata['morphology_mip']["scale0"].image.shape
    #h, w = sdata
    
    # Check if the image is already below the maximum number of pixels
    if h * w <= max_n_pixels:
        return None
    
    # Initialize with square crop
    h_crop = w_crop = int(np.sqrt(max_n_pixels))
    
    # Adjust crop if necessary to fit the image
    if h_crop > h:
        h_crop = h
        w_crop = int(max_n_pixels / h_crop)
    elif w_crop > w:
        w_crop = w
        h_crop = int(max_n_pixels / w_crop)
        
    # Center the crop
    h_offset = (h - h_crop) // 2
    w_offset = (w - w_crop) // 2
        
    crop = [[h_offset, h_offset + h_crop], [w_offset, w_offset + w_crop]]
        
    return crop

def rechunk_sdata(sdata, CHUNK_SIZE=1024):
    """Rechunk the sdata to the given chunk size
    
    Arguments
    ---------
    sdata: spatialdata.SpatialData
        The spatial data to rechunk
    CHUNK_SIZE: int
        The chunk size to rechunk to
        
    """
    
    for key in list(sdata.images.keys()):
        image = sdata.images[key]
        coords = list(image["scale0"].coords.keys())
        rechunk_strategy = {c: CHUNK_SIZE for c in coords}
        if "c" in coords:
            rechunk_strategy["c"] = image["scale0"]["image"].chunks[0][0]
        image = image.chunk(rechunk_strategy)
        sdata.images[key] = image
        
    for key in list(sdata.labels.keys()):
        label_image = sdata.labels[key]
        coords = list(label_image.coords.keys())
        rechunk_strategy = {c: CHUNK_SIZE for c in coords}
        label_image = label_image.chunk(rechunk_strategy)
        sdata.labels[key] = label_image


def subsample_adata_group_balanced(adata, group_key, n_samples, seed=0):
    """Subsample adata to a given number of samples, removing cells from large groups first
    
    Arguments
    ---------
    adata: anndata.AnnData
        The adata to subsample
    group_key: str
        The key in adata.obs to group by
    n_samples: int
        The number of samples to subsample to
    seed: int
        The seed to use for the random subsampling
        
    Returns
    -------
    pd.Series
        The series with the subsample information (boolean, True if the cell is in the subsample). 
        Series index is the same as adata.obs_names.
    """
    
    np.random.seed(seed)
    
    # Get the number of cells per group
    n_cells = adata.obs[group_key].value_counts().sort_values(ascending=True)
    
    if n_cells.sum() <= n_samples:
        all_obs_df = adata.obs.copy()
        all_obs_df["in_subsample"] = True
        return all_obs_df["in_subsample"]
    
    n_celltypes = len(n_cells)
    
    # Find out which groups to subsample from
    df = pd.DataFrame({"n_cells": n_cells, "sum": 0, "n_samples":0}, dtype=int)
    subsample_from_idx = n_celltypes
    tmp = np.zeros(n_celltypes, dtype=int)
    for i in range(n_celltypes):
        tmp[i] = df.iloc[:i]["n_cells"].sum()
        tmp[i] += (n_celltypes - i) * df.iloc[i]["n_cells"]
        if tmp[i] >= n_samples:
            subsample_from_idx = i
            break
    df["sum"] = tmp

    # Get number of samples per group
    n_samples_no_sampling = df.iloc[:subsample_from_idx]["n_cells"].sum()
    n_samples_to_subsample = n_samples - n_samples_no_sampling
    n_samples_per_group = n_samples_to_subsample // (n_celltypes - subsample_from_idx)
    n_samples_per_group_remainder = n_samples_to_subsample % (n_celltypes - subsample_from_idx)
    n_samples = np.zeros(n_celltypes, dtype=int)
    for i in range(subsample_from_idx):
        n_samples[i] = df.iloc[i]["n_cells"]
    for i in range(subsample_from_idx, n_celltypes):
        n_samples[i] = n_samples_per_group
        if n_samples_per_group_remainder > 0:
            n_samples[i] += 1
            n_samples_per_group_remainder -= 1
    df["n_samples"] = n_samples
    
    # Subsample from the selected groups
    mask_df = adata.obs[[group_key]].copy()
    mask_df["in_subsample"] = False
    for i in range(subsample_from_idx):
        ct = df.index[i]
        mask_df.loc[mask_df[group_key] == ct, "in_subsample"] = True
    for i in range(subsample_from_idx, n_celltypes):
        ct = df.index[i]
        ct_obs = mask_df.loc[mask_df[group_key] == ct].index
        ct_obs_subsample = np.random.choice(ct_obs, size=df.iloc[i]["n_samples"], replace=False)
        mask_df.loc[ct_obs_subsample, "in_subsample"] = True
        
    return mask_df["in_subsample"]



# Load the single-cell data
adata = ad.read_h5ad(par["input_sc"])

# Load the spatial data
sdata = sd.read_zarr(par["input_sp"])

# Subset single-cell data if it is too large
N_MAX_SC = 120000
if adata.n_obs > N_MAX_SC:
    adata = adata[subsample_adata_group_balanced(adata, "cell_type", N_MAX_SC, seed=0)]

# Subset single-cell and spatial data to shared genes
sp_genes = sdata['transcripts']['feature_name'].unique().compute().tolist()
sc_genes = adata.var["feature_name"].unique().tolist()
shared_genes = list(set(sp_genes) & set(sc_genes))
sdata['transcripts'] = sdata['transcripts'].loc[sdata['transcripts']['feature_name'].isin(shared_genes)]
adata = adata[:,adata.var["feature_name"].isin(shared_genes)].copy()

# Use feature names for adata instead of feature ids. convert to str
adata.var.reset_index(inplace=True, drop=True)
adata.var_names = adata.var["feature_name"].values.astype(str).tolist()

# store metadata to adata and sdata uns
metadata_uns_cols = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for col in metadata_uns_cols:
    orig_col = f"orig_{col}"
    if orig_col in adata.uns:
        adata.uns[orig_col] = adata.uns[col]
    adata.uns[col] = par[col]
    if orig_col in sdata.table.uns:
        sdata.table.uns[orig_col] = sdata.table.uns[col]
    sdata.table.uns[col] = par[col]

# Correct the feature_key attribute in sdata if needed
# NOTE: it would have been better to do this in the loader scripts, but this way the datasets don't need to be re-downloaded
if "feature_key" in sdata['transcripts'].attrs["spatialdata_attrs"]:
    feature_key = sdata['transcripts'].attrs["spatialdata_attrs"]["feature_key"]
    if feature_key != "feature_name":
        sdata['transcripts'].attrs["spatialdata_attrs"]["feature_key"] = "feature_name"

# Crop datasets that are too large
crop_coords = get_crop_coords(sdata)
if crop_coords is not None:
    sdata_output = sdata.query.bounding_box(
        axes=["y", "x"],
        min_coordinate=[crop_coords[0][0], crop_coords[1][0]],
        max_coordinate=[crop_coords[0][1], crop_coords[1][1]],
        target_coordinate_system="global",
        filter_table=True,
    )
    rechunk_sdata(sdata_output) #NOTE: rechunking currently needed (https://github.com/scverse/spatialdata/issues/929)
else:
    sdata_output = sdata
    
# Save the single-cell data
adata.write_h5ad(par["output_sc"], compression="gzip")

# remove directory if it exists
if os.path.exists(par["output_sp"]):
    shutil.rmtree(par["output_sp"])

# Save the spatial data
sdata_output.write(par["output_sp"], overwrite=True)
