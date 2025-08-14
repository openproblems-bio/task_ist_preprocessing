#!/usr/bin/env python3

import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import shutil
import os
import sys
import logging
from spatialdata.models import Labels2DModel
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

## VIASH START
par = {
    "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/dataset.zarr",
    "output": "segmentation.zarr",
    "transcripts_key": "transcripts",
    "shapes_key": "cell_boundaries",
    "images_key": "morphology_mip",
    "patch_width": 1200,
    "patch_overlap": 50,
    "transcript_patch_width": 200,
    "mean_cell_diameter": 15.0,
    "max_cell_radius": 25.0,
    "alpha": 0.5,
    "min_rna_per_cell": 5,
    "gene_column": "feature_name",
    "norm_vector": False,
    "allow_disconnected_polygon": True
}
## VIASH END

def validate_input_data(sdata):
    """Validate that required data components are present."""
    required_components = {
        'points': par['transcripts_key'],
        'images': par['images_key']
    }
    
    for component_type, key in required_components.items():
        if component_type not in sdata or key not in getattr(sdata, component_type):
            raise ValueError(f"Missing required {component_type} component: {key}")
    
    logger.info(f"Input data validated successfully")
    logger.info(f"Data components: {list(sdata)}")

def main():
    try:
        # Load input data
        logger.info(f"Loading data from {par['input']}")
        sdata_raw = sd.read_zarr(par["input"])
        
        # Validate input data
        validate_input_data(sdata_raw)
        
        # Check if prior segmentation exists
        prior_segmentation = "labels" in sdata_raw and par.get("segmentation_key") in sdata_raw.labels
        
        # Prepare SpatialData object for ComSeg
        components = {
            "points": {par["transcripts_key"]: sdata_raw.points[par["transcripts_key"]]},
            "images": {par["images_key"]: sdata_raw.images[par["images_key"]]}
        }
        
        # Add shapes if available
        if "shapes" in sdata_raw and par["shapes_key"] in sdata_raw.shapes:
            components["shapes"] = {par["shapes_key"]: sdata_raw.shapes[par["shapes_key"]]}
        
        # Add prior segmentation if available
        if prior_segmentation:
            components["labels"] = {"segmentation": sdata_raw.labels[par["segmentation_key"]]}
        
        # Create working SpatialData object
        sdata = sd.SpatialData(**components)
        
        logger.info(f"Prepared SpatialData for ComSeg processing")
        
        # Create image patches for analysis
        logger.info("Creating image patches")
        sopa.make_image_patches(
            sdata, 
            patch_width=par["patch_width"], 
            patch_overlap=par["patch_overlap"]
        )
        
        # Create transcript patches
        logger.info("Creating transcript patches")
        transcript_patch_args = {
            "sdata": sdata,
            "write_cells_centroids": True,
            "patch_width": par["transcript_patch_width"],
        }
        
        # Add prior shapes if available
        if par["shapes_key"] in components.get("shapes", {}):
            transcript_patch_args["prior_shapes_key"] = par["shapes_key"]
        
        sopa.make_transcript_patches(**transcript_patch_args)
        
        # Configure ComSeg parameters
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
        
        logger.info(f"Running ComSeg segmentation with config: {config}")
        
        # Run ComSeg segmentation
        sopa.segmentation.comseg(sdata, config)
        
        # Clean up temporary boundaries if created
        if "comseg_boundaries" in sdata.shapes:
            del sdata.shapes["comseg_boundaries"]
        
        # Create output SpatialData with segmentation results
        sd_output = sd.SpatialData()
        
        # Add segmentation labels if they exist
        if "segmentation" in sdata.labels:
            sd_output.labels["segmentation"] = sdata.labels["segmentation"]
        
        # Add transcript assignments if available
        if par["transcripts_key"] in sdata.points:
            transcripts = sdata.points[par["transcripts_key"]]
            if "cell_id" in transcripts.columns:
                # Create a minimal table with cell assignments
                unique_cells = transcripts["cell_id"].dropna().unique()
                if len(unique_cells) > 0:
                    # Create basic cell metadata
                    cell_df = pd.DataFrame({"cell_id": unique_cells})
                    cell_df.index = cell_df.index.astype(str)
                    
                    # Create minimal var data if available from original data
                    if "tables" in sdata_raw and "table" in sdata_raw.tables:
                        var_data = sdata_raw.tables["table"].var.copy()
                    else:
                        # Create minimal var data
                        genes = transcripts[par["gene_column"]].unique()
                        var_data = pd.DataFrame(index=genes)
                    
                    # Create AnnData table
                    adata = ad.AnnData(
                        obs=cell_df,
                        var=var_data
                    )
                    sd_output.tables["table"] = adata
        
        # Validate output
        if "segmentation" not in sd_output.labels:
            raise ValueError("ComSeg failed to produce segmentation results")
        
        logger.info("ComSeg segmentation completed successfully")
        logger.info(f"Output components: {list(sd_output)}")
        
        # Write output
        logger.info(f"Writing output to {par['output']}")
        if os.path.exists(par["output"]):
            shutil.rmtree(par["output"])
        
        sd_output.write(par["output"])
        logger.info("Output written successfully")
        
    except Exception as e:
        logger.error(f"ComSeg segmentation failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()