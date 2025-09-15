#!/usr/bin/env python3

import spatialdata as sd
import scanpy as sc
import numpy as np
import pandas as pd
import tifffile
import cv2
import natsort
import os 
import sys
import logging
import tempfile
from pathlib import Path
from bidcell import BIDCellModel

## VIASH START
par = {
    'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
    'output': 'output.zarr',
    'single_cell_ref': None,
    'max_overlaps_pos': 4,
    'max_overlaps_neg': 15,
    'model_epochs': 10,
    'min_cell_size': 15
}
## VIASH END

def generate_markers(ref_df, max_overlaps_pos=4, max_overlaps_neg=15):
    """Generate positive and negative marker genes for cell types."""
    n_genes = ref_df.shape[1] - 3  # Exclude ct_idx, cell_type, atlas columns
    cell_types = natsort.natsorted(list(set(ref_df["cell_type"].tolist())))
    n_cell_types = len(cell_types)

    ref_expr = ref_df.iloc[:, :n_genes].to_numpy()
    gene_names = ref_df.columns[:n_genes]
    ct_idx = ref_df["ct_idx"].to_numpy()

    # Generate negative markers (genes with low expression)
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

    # Generate positive markers (genes with high expression)
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

def main():
    logging.basicConfig(level=logging.INFO)
    
    # Load spatial data
    logging.info(f"Loading spatial data from {par['input']}")
    sdata = sd.read_zarr(par['input'])
    
    # Log data characteristics
    logging.info(f"Loaded spatial data with components: {list(sdata)}")
    
    # Get gene list from spatial transcripts
    sdata_genes = sdata['transcripts']["feature_name"].unique().compute().sort_values().tolist()
    logging.info(f"Found {len(sdata_genes)} unique genes in spatial data")
    
    # Create temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        try:
            # Extract DAPI image for BIDCell
            logging.info("Extracting morphology image")
            img = sdata["morphology_mip"]["scale0"]["image"].values
            img = np.squeeze(img)
            
            morphology_path = temp_path / "morphology_mip_pyramidal.tiff"
            with tifffile.TiffWriter(morphology_path, bigtiff=True) as tiff:
                tiff.write(img, photometric="minisblack", resolution=(1, 1))
            
            # Process single-cell reference data if provided
            # First check if there's an scRNA-seq reference in the spatial data itself
            if 'scrnaseq_reference' in sdata.tables:
                logging.info("Using scRNA-seq reference from input spatial data")
                adata = sdata.tables['scrnaseq_reference']
            elif par.get('single_cell_ref') and os.path.exists(par['single_cell_ref']):
                logging.info(f"Loading single-cell reference from {par['single_cell_ref']}")
                adata = sc.read_h5ad(par['single_cell_ref'])
            else:
                logging.info("No single-cell reference found, using scrnaseq_reference.h5ad from test data")
                # Try to use the scrnaseq_reference.h5ad that should be in the same directory as raw_ist.zarr
                input_dir = os.path.dirname(par['input'])
                ref_path = os.path.join(input_dir, 'scrnaseq_reference.h5ad')
                if os.path.exists(ref_path):
                    logging.info(f"Found scRNA-seq reference at {ref_path}")
                    adata = sc.read_h5ad(ref_path)
                else:
                    adata = None
                
            if adata is not None:
                # Filter to shared genes
                shared_genes = [g for g in sdata_genes if g in adata.var["feature_name"].values]
                logging.info(f"Found {len(shared_genes)} shared genes between spatial and scRNA-seq data")
                
                if len(shared_genes) == 0:
                    raise ValueError("No shared genes found between spatial and single-cell reference data")
                
                adata = adata[:, adata.var["feature_name"].isin(shared_genes)]
                adata.var_names = adata.var["feature_name"].astype(str)
                
                # Create reference dataframe for BIDCell
                # Use normalized layer if available, otherwise X
                if "normalized" in adata.layers:
                    expr_data = adata[:, shared_genes].layers["normalized"].toarray()
                else:
                    expr_data = adata[:, shared_genes].X.toarray()
                
                sc_ref = pd.DataFrame(
                    data=expr_data,
                    columns=shared_genes,
                    index=range(adata.n_obs)
                )
                
                # Add cell type information
                if 'cell_type' not in adata.obs.columns:
                    logging.warning("No 'cell_type' column found in reference data, using dummy cell type")
                    adata.obs['cell_type'] = 'Unknown'
                
                cell_type_col = adata.obs['cell_type'].astype('category')
                sc_ref["ct_idx"] = cell_type_col.cat.codes.values
                sc_ref["cell_type"] = cell_type_col.values
                sc_ref["atlas"] = "custom"
                
                # Save reference data
                scref_path = temp_path / "scref.csv"
                sc_ref.to_csv(scref_path)
                
                # Generate marker files
                logging.info("Generating positive and negative marker genes")
                df_pos, df_neg = generate_markers(
                    sc_ref, 
                    max_overlaps_pos=par['max_overlaps_pos'],
                    max_overlaps_neg=par['max_overlaps_neg']
                )
                
                pos_marker_path = temp_path / "pos_marker.csv"
                neg_marker_path = temp_path / "neg_marker.csv"
                df_pos.to_csv(pos_marker_path)
                df_neg.to_csv(neg_marker_path)
                
                # Filter transcripts to shared genes
                transcript = sdata["transcripts"].compute()
                transcript_filtered = transcript[transcript["feature_name"].isin(shared_genes)]
                
            else:
                logging.warning("No single-cell reference provided, using all genes")
                transcript_filtered = sdata["transcripts"].compute()
                shared_genes = sdata_genes
            
            # Save transcript data for BIDCell
            transcript_path = temp_path / "transcript.csv.gz"
            pd.DataFrame(transcript_filtered).to_csv(transcript_path, compression='gzip')
            
            # Create BIDCell configuration file
            config = {
                'data_path': str(temp_path),
                'morphology_path': str(morphology_path),
                'transcript_path': str(transcript_path),
                'epochs': par['model_epochs'],
                'min_cell_size': par['min_cell_size']
            }
            
            if adata is not None:
                config.update({
                    'scref_path': str(scref_path),
                    'pos_marker_path': str(pos_marker_path),
                    'neg_marker_path': str(neg_marker_path)
                })
            
            config_path = temp_path / "bidcell_config.yaml"
            import yaml
            with open(config_path, 'w') as f:
                yaml.dump(config, f)
            
            # Run BIDCell
            logging.info("Running BIDCell segmentation")
            model = BIDCellModel(str(config_path))
            model.run_pipeline()
            
            # Process BIDCell output
            logging.info("Processing BIDCell output")
            dapi_image = tifffile.imread(morphology_path)
            
            # Look for BIDCell output file (adjust name based on actual output)
            output_files = list(temp_path.glob("*connected.tif"))
            if not output_files:
                output_files = list(temp_path.glob("segmentation*.tif"))
            
            if not output_files:
                raise FileNotFoundError("BIDCell segmentation output not found")
            
            segmentation_mask = tifffile.imread(output_files[0])
            h_dapi, w_dapi = dapi_image.shape
            
            # Resize segmentation to match DAPI image
            segmentation_resized = cv2.resize(
                segmentation_mask.astype('float32'), 
                (w_dapi, h_dapi), 
                interpolation=cv2.INTER_NEAREST
            )
            segmentation_resized = segmentation_resized.astype(np.uint32)
            
            # Create output SpatialData
            logging.info("Creating output SpatialData")
            
            # Prepare images
            image_with_channel = np.expand_dims(dapi_image, axis=0)
            images = sd.models.Image2DModel.parse(image_with_channel, dims=('c', 'y', 'x'))
            
            # Prepare labels (segmentation)
            labels = sd.models.Labels2DModel.parse(segmentation_resized, dims=('y', 'x'))
            
            # Prepare points (transcripts)
            transcript_df = pd.DataFrame(transcript_filtered)
            transcript_df['x'] = transcript_df['x'].astype(float)
            transcript_df['y'] = transcript_df['y'].astype(float)
            transcript_df['z'] = transcript_df['z'].astype(float)
            transcript_df['feature_name'] = transcript_df['feature_name'].astype(str)
            points = sd.models.PointsModel.parse(transcript_df)
            
            # Create output SpatialData object
            output_sdata = sd.SpatialData(
                images={'morphology_mip': images},
                labels={'segmentation': labels},
                points={'transcripts': points}
            )
            
            # Write output
            logging.info(f"Writing output to {par['output']}")
            output_sdata.write(par['output'], overwrite=True)
            
            logging.info("BIDCell segmentation completed successfully")
            
        except Exception as e:
            logging.error(f"BIDCell segmentation failed: {str(e)}")
            sys.exit(1)

if __name__ == "__main__":
    main()
