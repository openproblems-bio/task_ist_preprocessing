#!/bin/bash

set -e

SP_DIR=resources_test/common/2023_10x_mouse_brain_xenium_rep1
SC_DIR=resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2
OUT_DIR="resources_test/task_ist_preprocessing/mouse_brain_combined"

rm -rf $OUT_DIR
mkdir -p $OUT_DIR


# run dataset preprocessor
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input_sc $SC_DIR/dataset.h5ad \
  --input_sp $SP_DIR/dataset.zarr \
  --output_sc $OUT_DIR/scrnaseq_reference.h5ad \
  --output_sp $OUT_DIR/raw_ist.zarr

# run a segmentation method
viash run src/methods_segmentation/custom_segmentation/config.vsh.yaml -- \
  --input $OUT_DIR/raw_ist.zarr \
  --labels_key cell_labels \
  --output $OUT_DIR/segmentation.zarr

# run an assignment method
viash run src/methods_transcript_assignment/basic_transcript_assignment/config.vsh.yaml -- \
  --input_ist $OUT_DIR/raw_ist.zarr \
  --input_segmentation $OUT_DIR/segmentation.zarr \
  --transcripts_key transcripts \
  --coordinate_system global \
  --output $OUT_DIR/transcript_assignments.zarr

# run a count aggregation method
viash run src/methods_count_aggregation/basic_count_aggregation/config.vsh.yaml -- \
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
  --input_transcript_assignments $OUT_DIR/transcript_assignments.zarr \
  --input_scrnaseq_reference $OUT_DIR/scrnaseq_reference.h5ad \
  --output $OUT_DIR/spatial_with_cell_types.h5ad

# run a gene efficiency correction method
viash run src/methods_expression_correction/gene_efficiency_correction/config.vsh.yaml -- \
  --input_spatial_with_cell_types $OUT_DIR/spatial_with_cell_types.h5ad \
  --input_scrnaseq_reference $OUT_DIR/scrnaseq_reference.h5ad \
  --output $OUT_DIR/spatial_corrected_counts.h5ad

# run a QC filter method
viash run src/methods_qc_filter/basic_qc_filter/config.vsh.yaml -- \
  --input $OUT_DIR/spatial_corrected_counts.h5ad \
  --output $OUT_DIR/spatial_qc_col.h5ad

# run a metric
viash run src/metrics/similarity/config.vsh.yaml -- \
  --input $OUT_DIR/spatial_corrected_counts.h5ad \
  --input_qc_col $OUT_DIR/spatial_qc_col.h5ad \
  --input_sc $OUT_DIR/scrnaseq_reference.h5ad \
  --input_transcript_assignments $OUT_DIR/transcript_assignments.zarr \
  --output $OUT_DIR/score.h5ad

# create a state file
cat >> $OUT_DIR/state.yaml <<EOL
id: mouse_brain_combined
output_sp: !file raw_ist.zarr
output_sc: !file scrnaseq_reference.h5ad
output_segmentation: !file segmentation.zarr
output_transcript_assignments: !file transcript_assignments.zarr
output_spatial_aggregated_counts: !file spatial_aggregated_counts.h5ad
output_cell_volumes: !file cell_volumes.h5ad
output_spatial_normalized_counts: !file spatial_normalized_counts.h5ad
output_spatial_with_cell_types: !file spatial_with_cell_types.h5ad
output_spatial_corrected_counts: !file spatial_corrected_counts.h5ad
output_spatial_qc_col: !file spatial_qc_col.h5ad
output_score: !file score.h5ad
EOL

# sync test resources
aws s3 sync --profile op \
  "resources_test/task_ist_preprocessing/mouse_brain_combined" \
  "s3://openproblems-data/resources_test/task_ist_preprocessing/mouse_brain_combined" \
  --delete --dryrun
