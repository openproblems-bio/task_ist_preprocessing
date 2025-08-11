import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import shutil

# Testing
par = {
    "input": "/home/anbui/cmscb/testdata/common-20250310T113316Z-001/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr",
    "output": "output.zarr",
    "transcripts_key": "transcripts",
    "segmentation_key": "cell_labels",  # Change if needed
    "shapes_key": "cell_boundaries",  # Required for sopa
}

sdata_raw = sd.read_zarr(par["input"])

print(sdata_raw)

print(sdata_raw.labels)

prior_segmentation = par.get("segmentation_key") in sdata_raw.labels

if prior_segmentation:
    # Read segmentation labels
    sdata_segm = sdata_raw.labels[par["segmentation_key"]]

    # Create a new SpatialData object with segmentation and transcripts
    sdata = sd.SpatialData(
        labels={"segmentation": sdata_segm},
        points={"transcripts": sdata_raw.points[par["transcripts_key"]]},
        shapes={"cell_boundaries": sdata_raw["cell_boundaries"]},  # Added cell boundaries
        images={"morphology_mip": sdata_raw["morphology_mip"]},
       
    )
else:
    # Create a new SpatialData object with only transcripts
    sdata = sd.SpatialData(
        points={"transcripts": sdata_raw.points[par["transcripts_key"]]},
        shapes={"cell_boundaries": sdata_raw["cell_boundaries"]},  # Added cell boundaries
        images={"morphology_mip": sdata_raw["morphology_mip"]},
       
    )

    sdata.write(par['output'])

    # Create patches for analysis
sopa.make_image_patches(sdata, patch_width=1200, patch_overlap=50)

# Make transcript patches with your dataset's structure
sopa.make_transcript_patches(
    sdata,
    prior_shapes_key="cell_boundaries",
    write_cells_centroids=True,
    patch_width=200,
)

# Configure comseg with dataset's structure
config = {
    "dict_scale": {"x": 1, "y": 1, "z": 1},  # spot coordinates already in µm
    "mean_cell_diameter": 15,
    "max_cell_radius": 25,
    "norm_vector": False,
    "alpha": 0.5, 
    "allow_disconnected_polygon": True,
    "min_rna_per_cell": 5,  
    "gene_column": "feature_name",  
}

# Run comseg segmentation
sopa.segmentation.comseg(sdata, config)

del sdata["comseg_boundaries"]

print(sdata["transcripts"].columns)

cell_id_col = sdata["transcripts"]["cell_id"]

sdata.tables["table"]=ad.AnnData(obs=pd.DataFrame({"cell_id":cell_id_col}), var=sdata_raw.tables["table"].var[[]])


temp_output = par["output"] + "_temp"

# Write to a temporary location
sdata.write(temp_output, overwrite=True)

# Remove the old output and rename
shutil.rmtree(par["output"], ignore_errors=True)  # Delete old output
shutil.move(temp_output, par["output"])  # Rename temp output to original

sdata_new = sd.SpatialData(
    points=sdata.points,  
    tables=sdata.tables   
)

output_path = "cleaned_output.zarr" 
sdata_new.write(output_path, overwrite=True)

print(sdata_new)

