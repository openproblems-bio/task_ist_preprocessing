import spatialdata as sd
import spatialdata_plot as pl
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import cv2
import dask.dataframe as dd
import scanpy as sc
import pandas as pd
import natsort
import os 
from bidcell import BIDCellModel

## VIASH START
par = {
    'input': "../task_ist_preprocessing/resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
    'temp': './temp/bidcell/',
    'output': 'output.zarr',
    'single_cell_ref': None,
    'max_overlaps_pos': 4,
    'max_overlaps_neg': 15,
    'model_epochs': 10,
    'min_cell_size': 15
}
## VIASH END

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

if not os.path.exists(par['temp']):
    os.makedirs(par['temp'])


sdata = sd.read_zarr(par['input'])
sdata_genes = sdata['transcripts']["feature_name"].unique().compute().sort_values().tolist()


image_pyramid = []
img = sdata["morphology_mip"]["/scale0"]["image"].values  # Convert dask array to numpy
img = np.squeeze(img)  # Remove singleton channel dimension (c:1)
image_pyramid.append(img)



# Save the TIFF file in the temporary directory
with tifffile.TiffWriter(f"{par['temp']}morphology_mip_pyramidal.tiff", bigtiff=True) as tiff:
    for img in image_pyramid:
        tiff.write(img, photometric="minisblack", resolution=(1, 1))



adata = sc.read_h5ad(par['input_scrnaseq_reference'])
shared_genes = [g for g in sdata_genes if g in adata.var["feature_name"].values]
adata = adata[:, adata.var["feature_name"].isin(shared_genes)]
adata.var_names = adata.var["feature_name"].astype(str)


# Make scref.csv
sc_ref = pd.DataFrame(
    data=adata.layers["normalized"].toarray(),
    columns=shared_genes,
    index=range(adata.n_obs)
)


celltypes = adata.obs['cell_type'].unique().tolist()
cell_type_col = adata.obs['cell_type'].astype('category')
sc_ref["ct_idx"] = cell_type_col.cat.codes.values
sc_ref["cell_type"] = cell_type_col.values
sc_ref["atlas"] = "custom"
sc_ref.to_csv(f"{par['temp']}scref.csv")


transcript = sdata["transcripts"].compute()
transcript = pd.DataFrame(transcript)
transcript[transcript["feature_name"].isin(shared_genes)].to_csv(
    f"{par['temp']}transcript.csv.gz", compression="gzip"
)




# generate positive and negative marker files 
df_pos, df_neg = generate_markers(sc_ref, max_overlaps_pos=4, max_overlaps_neg=15)
df_pos.to_csv(f"{par['temp']}/pos_marker.csv")
df_neg.to_csv(f"{par['temp']}/neg_marker.csv")


import yaml

config = {
    "cpus": 8, 
    "files": {
        "data_dir": par['temp'],
        "fp_dapi": f"{par['temp']}morphology_mip_pyramidal.tiff",
        "fp_transcripts": f"{par['temp']}transcript.csv.gz",
        "fp_ref": f"{par['temp']}scref.csv",
        "fp_pos_markers": f"{par['temp']}pos_marker.csv",
        "fp_neg_markers": f"{par['temp']}neg_marker.csv",
    },
    "nuclei_fovs": {
        "stitch_nuclei_fovs": False,
    },
    "nuclei": {
        "diameter": None,  # leave as None to automatically compute
    },
    "transcripts": {
        "shift_to_origin": True,
        "x_col": "x",
        "y_col": "y",
        "gene_col": "feature_name",
        "transcripts_to_filter": [
            "NegControlProbe_",
            "antisense_",
            "NegControlCodeword_",
            "BLANK_",
            "Blank-",
            "NegPrb"
        ],

    },
    "affine": {
        "target_pix_um": 1.0,
        "base_pix_x": 0.2125,
        "base_pix_y": 0.2125,
        "base_ts_x": 0.1,
        "base_ts_y": 0.1,
        "global_shift_x": 0,
        "global_shift_y": 0,
    },
    "model_params": {
        "name": "custom",
        "patch_size": 48,
        "elongated": list(sc_ref["cell_type"]),
    },
    "training_params": {
        "total_epochs": 1,
        "total_steps": 60,
        "ne_weight": 1.0,
        "os_weight": 1.0,
        "cc_weight": 1.0,
        "ov_weight": 1.0,
        "pos_weight": 1.0,
        "neg_weight": 1.0,
    },
    "testing_params": {
        "test_epoch": 1,
        "test_step": 60,
    },
    "experiment_dirs": {
        "dir_id": "last",
    },
}

# Save YAML file
with open(f"{par['temp']}testdata.yaml", "w") as f:
    yaml.dump(config, f, sort_keys=False)




#for i in sc_ref.columns:
#    print(i, flush = True)
#    # Setting up and running BIDCell
model = BIDCellModel(f"{par['temp']}testdata.yaml")
model.run_pipeline() 




    # Analysis and visualisation of BIDcell output
#dapi_image = tifffile.imread("morphology_mip_pyramidal.tiff")
#segmentation_mask = tifffile.imread("epoch_10_step_60_connected.tif")
#h_dapi, w_dapi = dapi_image.shape

#segmentation_mask_resized = cv2.resize(segmentation_mask.astype('float32'), (w_dapi, h_dapi), interpolation=cv2.INTER_NEAREST)
#segmentation_mask_resized = segmentation_mask_resized.astype(np.uint32)
#segmentation_mask_resized = segmentation_mask_resized.transpose(1, 0)
#tifffile.imwrite("bidcellresult_resized.tif", segmentation_mask_resized)

    # creating bidcelloutput.zarr
#image = tifffile.imread("morphology_mip_pyramidal.tiff")
#image_with_channel = np.expand_dims(image, axis=0)
#label_image = tifffile.imread("bidcellresult_resized.tif")
#labels = sd.models.Labels2DModel.parse(label_image, dims=('y', 'x'))
    
#transcript_processed = pd.read_csv("data/transcript.csv.gz")
#transcript_processed['x'] = transcript_processed['x'].astype(float)
#transcript_processed['y'] = transcript_processed['y'].astype(float)
#transcript_processed['z'] = transcript_processed['z'].astype(float)
#transcript_processed['feature_name'] = transcript_processed['feature_name'].astype(str)

#images = sd.models.Image2DModel.parse(image_with_channel, dims=('c', 'x', 'y'))
#labels = sd.models.Labels2DModel.parse(label_image, dims=('y', 'x'))
#points = sd.models.PointsModel.parse(transcript_processed)

outputsdata = sd.SpatialData()
#        images={'DAPI': images},
#        labels={'segmentation_mask_labels': labels},  
#        points={'transcripts': points}  
#    )
outputsdata.write(par['output'], overwrite=True)