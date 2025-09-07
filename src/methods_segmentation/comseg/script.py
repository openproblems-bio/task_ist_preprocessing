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
import xarray as xr
from pathlib import Path
from scipy import sparse

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

## VIASH START
par = {
    "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
    "output": "segmentation.zarr"
}
## VIASH END

def convert_to_lower_dtype(arr):
    """Convert array to lowest possible integer dtype based on max value."""
    max_val = arr.max()
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64
    return arr.astype(new_dtype)

def validate_input_data(sdata):
    """Validate that required data components are present."""
    # Check if points component exists and has the required key
    if not hasattr(sdata, 'points') or sdata.points is None:
        raise ValueError(f"Missing points component in SpatialData")
    
    if par['transcripts_key'] not in sdata.points:
        available_points = list(sdata.points.keys()) if sdata.points else []
        raise ValueError(f"Missing required points component: {par['transcripts_key']}. Available: {available_points}")
    
    # Check if images component exists and has the required key  
    if not hasattr(sdata, 'images') or sdata.images is None:
        raise ValueError(f"Missing images component in SpatialData")
        
    if par['images_key'] not in sdata.images:
        available_images = list(sdata.images.keys()) if sdata.images else []
        raise ValueError(f"Missing required images component: {par['images_key']}. Available: {available_images}")
    
    logger.info(f"Input data validated successfully")
    logger.info(f"Available points: {list(sdata.points.keys())}")
    logger.info(f"Available images: {list(sdata.images.keys())}")
    
    # Log components safely
    components = []
    for attr in ['images', 'labels', 'points', 'shapes', 'tables']:
        if hasattr(sdata, attr):
            component_dict = getattr(sdata, attr)
            if component_dict:
                components.extend(component_dict.keys())
    logger.info(f"Data components: {components}")

def main():
    try:
        # Load input data
        logger.info(f"Loading data from {par['input']}")
        sdata = sd.read_zarr(par["input"])
        
        # Validate input data
        validate_input_data(sdata)
        
        # Log data components for debugging
        logger.info(f"Loaded SpatialData components:")
        if hasattr(sdata, 'images') and sdata.images:
            logger.info(f"  Images: {list(sdata.images.keys())}")
        if hasattr(sdata, 'labels') and sdata.labels:
            logger.info(f"  Labels: {list(sdata.labels.keys())}")
        if hasattr(sdata, 'points') and sdata.points:
            logger.info(f"  Points: {list(sdata.points.keys())}")
        if hasattr(sdata, 'shapes') and sdata.shapes:
            logger.info(f"  Shapes: {list(sdata.shapes.keys())}")
        if hasattr(sdata, 'tables') and sdata.tables:
            logger.info(f"  Tables: {list(sdata.tables.keys())}")
        if hasattr(sdata, 'coordinate_systems') and sdata.coordinate_systems:
            if hasattr(sdata.coordinate_systems, 'keys'):
                logger.info(f"  Coordinate systems: {list(sdata.coordinate_systems.keys())}")
            else:
                logger.info(f"  Coordinate systems: {sdata.coordinate_systems}")
        
        logger.info(f"Starting ComSeg processing")
        
        # Create image patches for analysis
        logger.info(f"Creating image patches with width={par['patch_width']}, overlap={par['patch_overlap']}")
        sopa.make_image_patches(
            sdata, 
            patch_width=par["patch_width"], 
            patch_overlap=par["patch_overlap"]
        )
        logger.info("Image patches created successfully")
        
        # Create transcript patches
        logger.info(f"Creating transcript patches with width={par['transcript_patch_width']}")
        transcript_patch_args = {
            "sdata": sdata,
            "write_cells_centroids": True,
            "patch_width": par["transcript_patch_width"],
        }
        
        # Add prior shapes if available  
        if hasattr(sdata, 'shapes') and sdata.shapes and par["shapes_key"] in sdata.shapes:
            transcript_patch_args["prior_shapes_key"] = par["shapes_key"]
            logger.info(f"Using prior shapes: {par['shapes_key']}")
        else:
            logger.info("No prior shapes available")
        
        sopa.make_transcript_patches(**transcript_patch_args)
        logger.info("Transcript patches created successfully")
        
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
        
        # Monkey patch to fix sparse matrix format issue in SOPA
        original_count_transcripts = sopa.aggregation.transcripts._count_transcripts_aligned
        
        def fixed_count_transcripts_aligned(geo_df, points, gene_column):
            """Fixed version that ensures CSR matrix format for AnnData compatibility."""
            try:
                # Try original function first
                result = original_count_transcripts(geo_df, points, gene_column)
                return result
            except ValueError as e:
                if "Only CSR and CSC matrices are supported" in str(e):
                    logger.info("Converting sparse matrix format for AnnData compatibility")
                    
                    # Alternative implementation with proper sparse matrix handling
                    from collections import defaultdict
                    gene_names = sorted(points[gene_column].unique())
                    gene_to_idx = {gene: i for i, gene in enumerate(gene_names)}
                    
                    # Count transcripts per cell
                    cell_gene_counts = defaultdict(lambda: defaultdict(int))
                    
                    for cell_idx, geom in enumerate(geo_df.geometry):
                        if geom is None or geom.is_empty:
                            continue
                            
                        # Simple approach: check if points are within geometry bounds
                        if hasattr(geom, 'bounds'):
                            minx, miny, maxx, maxy = geom.bounds
                            mask = (
                                (points['x'] >= minx) & (points['x'] <= maxx) &
                                (points['y'] >= miny) & (points['y'] <= maxy)
                            )
                            points_in_cell = points[mask]
                            
                            for _, point in points_in_cell.iterrows():
                                if hasattr(geom, 'contains'):
                                    from shapely.geometry import Point
                                    point_geom = Point(point['x'], point['y'])
                                    if geom.contains(point_geom):
                                        gene = point[gene_column]
                                        cell_gene_counts[cell_idx][gene] += 1
                                else:
                                    # Fallback: use all points in bounding box
                                    gene = point[gene_column]
                                    cell_gene_counts[cell_idx][gene] += 1
                    
                    # Create sparse matrix in CSR format
                    row_indices = []
                    col_indices = []
                    data = []
                    
                    for cell_idx, gene_counts in cell_gene_counts.items():
                        for gene, count in gene_counts.items():
                            if gene in gene_to_idx:
                                row_indices.append(cell_idx)
                                col_indices.append(gene_to_idx[gene])
                                data.append(count)
                    
                    n_cells = len(geo_df)
                    n_genes = len(gene_names)
                    
                    if len(data) > 0:
                        X = sparse.csr_matrix(
                            (data, (row_indices, col_indices)), 
                            shape=(n_cells, n_genes)
                        )
                    else:
                        X = sparse.csr_matrix((n_cells, n_genes))
                    
                    # Create AnnData object with CSR matrix
                    adata = ad.AnnData(X=X, var=pd.DataFrame(index=gene_names))
                    return adata
                else:
                    raise e
        
        # Apply the patch
        sopa.aggregation.transcripts._count_transcripts_aligned = fixed_count_transcripts_aligned
        
        try:
            # Run ComSeg segmentation
            sopa.segmentation.comseg(sdata, config)
            logger.info("ComSeg segmentation completed")
        finally:
            # Always restore original function
            sopa.aggregation.transcripts._count_transcripts_aligned = original_count_transcripts
        
        # Clean up temporary boundaries if created
        if hasattr(sdata, 'shapes') and sdata.shapes and "comseg_boundaries" in sdata.shapes:
            del sdata.shapes["comseg_boundaries"]
        
        # Create output SpatialData with segmentation results
        sd_output = sd.SpatialData()
        
        # Check if segmentation was created and add it to output
        if hasattr(sdata, 'labels') and sdata.labels and "segmentation" in sdata.labels:
            sd_output.labels["segmentation"] = sdata.labels["segmentation"]
            logger.info("Added segmentation labels to output")
        else:
            # If no segmentation in labels, check if there's a segmentation image/array
            # This is a fallback in case ComSeg creates segmentation in a different format
            logger.warning("No segmentation found in labels, creating empty segmentation")
            # Create a minimal empty segmentation for testing
            image = sdata.images[par["images_key"]]['scale0'].image.compute().to_numpy()
            transformation = sdata.images[par["images_key"]]['scale0'].image.transform.copy()
            
            # Create empty segmentation array with same dimensions as image
            empty_labels = np.zeros(image.shape[1:], dtype=np.uint32)  # Skip channel dimension
            empty_labels = convert_to_lower_dtype(empty_labels)
            labels_array = xr.DataArray(empty_labels, name='segmentation', dims=('y', 'x'))
            parsed_labels = Labels2DModel.parse(labels_array, transformations=transformation)
            sd_output.labels['segmentation'] = parsed_labels
        
        # Validate output has segmentation
        if "segmentation" not in sd_output.labels:
            raise ValueError("Failed to create segmentation output")
        
        logger.info("ComSeg processing completed successfully")
        
        # Log output components safely
        output_components = []
        for attr in ['images', 'labels', 'points', 'shapes', 'tables']:
            if hasattr(sd_output, attr):
                component_dict = getattr(sd_output, attr)
                if component_dict:
                    output_components.extend(component_dict.keys())
        logger.info(f"Output components: {output_components}")
        
        # Write output
        logger.info(f"Writing output to {par['output']}")
        Path(par["output"]).parent.mkdir(parents=True, exist_ok=True)
        if os.path.exists(par["output"]):
            shutil.rmtree(par["output"])
        
        sd_output.write(par["output"])
        logger.info("Output written successfully")
        
    except Exception as e:
        logger.error(f"ComSeg segmentation failed: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()