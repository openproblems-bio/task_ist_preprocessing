import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import numpy as np

## VIASH START
par = {
    "input_ist": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
    "input_segmentation": "resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr",
    "transcripts_key": "transcripts",
    "coordinate_system": "global",
    "output": "temp/comseg/transcripts.zarr",
    
    "patch_width": 2000,
    "patch_overlap": 50,
    "transcript_patch_width": 200,
    "mean_cell_diameter": 15.0,
    "max_cell_radius": 25.0,
    "alpha": 0.5,
    "min_rna_per_cell": 5,
    "gene_column": "feature_name",
    "norm_vector": False,
    "allow_disconnected_polygon": True,
}
## VIASH END


# Read input files
print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])

# Convert the prior segmentation to polygons
sdata["segmentation_boundaries"] = sd.to_polygons(sdata_segm["segmentation"])
del sdata["segmentation_boundaries"]["label"] # make_transcript_patches will create a new label column and fails if one exists.

# Make patches
sopa.make_image_patches(sdata, image_key="morphology_mip", patch_width=par["patch_width"], patch_overlap=par["patch_overlap"])

transcript_patch_args = {
    "sdata": sdata,
    "write_cells_centroids": True,
    "patch_width": par["transcript_patch_width"],
    "prior_shapes_key": "segmentation_boundaries",
    "points_key": par["transcripts_key"],
}

sopa.make_transcript_patches(**transcript_patch_args)

# Run ComSeg
config = {
    "dict_scale": {"x": 1, "y": 1, "z": 1},
    "mean_cell_diameter": par["mean_cell_diameter"],
    "max_cell_radius": par["max_cell_radius"],
    "norm_vector": par["norm_vector"],
    "alpha": par["alpha"], 
    "allow_disconnected_polygon": par["allow_disconnected_polygon"],
    "min_rna_per_cell": par["min_rna_per_cell"],
    "gene_column": par["gene_column"],
}


# sopa.settings.parallelization_backend = 'dask' #TODO: get parallelization running.
sopa.segmentation.comseg(sdata, config)

# Assign transcripts to cell ids
sopa.spatial.assign_transcript_to_cell(
    sdata,
    points_key="transcripts",
    shapes_key="comseg_boundaries",
    key_added="cell_id",
    unassigned_value=0
)

# Create output SpatialData 

# Create objects for cells table
print('Creating objects for cells table', flush=True)
unique_cells = np.unique(sdata["transcripts"]["cell_id"])
zero_idx = np.where(unique_cells == 0)
if len(zero_idx[0]): 
    unique_cells=np.delete(unique_cells, zero_idx[0][0])

cell_id_col = pd.Series(unique_cells, name='cell_id', index=unique_cells)

# Create transcripts only sdata
print('Subsetting to transcripts cell id data', flush=True)
sdata_transcripts_only = sd.SpatialData(
    points={
        "transcripts": sdata['transcripts']
    },
    tables={
        "table": ad.AnnData(
          obs=pd.DataFrame(cell_id_col),
        )
    }
)


output_path = par['output']
sdata_transcripts_only.write(output_path, overwrite=True)


