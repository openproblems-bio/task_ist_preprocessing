#!/bin/bash

set -e

OUT_DIR="resources_test/task_ist_preprocessing/mouse_brain_combined"

# run dataset preprocessor
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input_scrnaseq resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad \
  --input_ist resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr \
  --output_scrnaseq $OUT_DIR/scrnaseq_reference.h5ad \
  --output_ist $OUT_DIR/raw_ist.zarr

# run a segmentation method
viash run src/methods_segmentation/custom/config.vsh.yaml -- \
  --input $OUT_DIR/raw_ist.zarr \
  --labels_key cell_labels \
  --output $OUT_DIR/segmentation.zarr

# run an assignment method
viash run src/methods_transcript_assignment/basic/config.vsh.yaml -- \
  --input_ist $OUT_DIR/raw_ist.zarr \
  --input_segmentation $OUT_DIR/segmentation.zarr \
  --transcripts_key transcripts \
  --coordinate_system global \
  --output $OUT_DIR/transcript_assignments.zarr

# run a count aggregation method
viash run src/methods_count_aggregation/basic/config.vsh.yaml -- \
  --input $OUT_DIR/transcript_assignments.zarr \
  --output $OUT_DIR/spatial_aggregated_counts.h5ad

# run a cell volume method
viash run src/methods_calculate_cell_volume/alpha_shapes/config.vsh.yaml -- \
  --input $OUT_DIR/transcript_assignments.zarr \
  --output $OUT_DIR/cell_volumes.h5ad

# run a normalization method
viash run src/methods_normalization/normalize_by_volume/config.vsh.yaml -- \
  --input_spatial_aggregated_counts $OUT_DIR/spatial_aggregated_counts.h5ad \
  --input_cell_volumes $OUT_DIR/cell_volumes.h5ad \
  --output $OUT_DIR/spatial_normalized_counts.h5ad

# run a cell type annotation method
viash run src/methods_cell_type_annotation/ssam/config.vsh.yaml -- \
  --input_spatial_normalized_counts $OUT_DIR/spatial_normalized_counts.h5ad \
  --input_transcripts $OUT_DIR/transcript_assignments.zarr \
  --input_scrnaseq_reference $OUT_DIR/scrnaseq_reference.h5ad \
  --output $OUT_DIR/spatial_with_cell_types.h5ad

# run a gene efficiency correction method
viash run src/methods_expression_correction/gene_efficiency_correction/config.vsh.yaml -- \
  --input_spatial_with_cell_types $OUT_DIR/spatial_with_cell_types.h5ad \
  --input_scrnaseq_reference $OUT_DIR/scrnaseq_reference.h5ad \
  --output $OUT_DIR/spatial_corrected_counts.h5ad