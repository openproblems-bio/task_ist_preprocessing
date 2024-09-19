#!/bin/bash

SP_DIR=resources_test/common/2023_10x_mouse_brain_xenium_rep1
SC_DIR=resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2
OUT_DIR="resources_test/task_ist_preprocessing/mouse_brain_combined"

viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input_sc resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad \
  --input_sp resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr \
  --output_sc $OUT_DIR/sc_normalised.h5ad \
  --output_sp $OUT_DIR/spatial_processed.zarr

viash run src/segmentation_methods/custom/config.vsh.yaml -- \
  --input $OUT_DIR/spatial_processed.zarr \
  --labels_key rep1_cell \
  --output $OUT_DIR/label_image.zarr

viash run src/assignment_methods/basic/config.vsh.yaml -- \
  --input $OUT_DIR/spatial_processed.zarr \
  --transcripts_key rep1_transcripts \
  --segmentation_input $OUT_DIR/label_image.zarr \
  --coordinate_system rep1_global \
  --output $OUT_DIR/assigned_transcripts.zarr

viash run src/cell_volume_methods/alpha_shapes/config.vsh.yaml -- \
  --input $OUT_DIR/assigned_transcripts.zarr \
  --output $OUT_DIR/cell_volumes.h5ad

viash run src/count_aggregation/basic/config.vsh.yaml -- \
  --input $OUT_DIR/assigned_transcripts.zarr \
  --output $OUT_DIR/raw_counts.h5ad

viash run src/normalisation_methods/normalise_by_volume/config.vsh.yaml -- \
  --input $OUT_DIR/raw_counts.h5ad \
  --input_volume $OUT_DIR/cell_volumes.h5ad \
  --output $OUT_DIR/norm_counts.h5ad

viash run src/celltype_annotation_methods/ssam/config.vsh.yaml -- \ 
    --input $OUT_DIR/norm_counts.h5ad \ 
    --input_transcripts $OUT_DIR/assigned_transcripts.zarr \ 
    --input_sc $OUT_DIR/sc_normalised.h5ad \ 
    --output $OUT_DIR/spatial_with_celltypes.h5ad

viash run src/expr_correction_methods/gene_efficiency_correction/config.vsh.yaml -- \ 
    --input $OUT_DIR/spatial_with_celltypes.h5ad \ 
    --input_sc $OUT_DIRsc_normalised.h5ad \ 
    --output $OUT_DIR/spatial_corrected.h5ad

# target/executable/segmentation_methods/custom/custom ---engine native --input resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr --labels_key rep1_cell --output $OUT_DIR/label_image.zarr
# target/executable/assignment_methods/basic/basic ---engine native --input resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr --transcripts_key rep1_transcripts --segmentation_input $OUT_DIR/label_image.zarr --coordinate_system rep1_global --output $OUT_DIR/assigned_transcripts.zarr
# target/executable/cell_volume_methods/alpha_shapes/alpha_shapes ---engine native --input $OUT_DIR/assigned_transcripts.zarr --output $OUT_DIR/cell_volumes.h5ad
# target/executable/count_aggregation/basic/basic ---engine native --input $OUT_DIR/assigned_transcripts.zarr --output $OUT_DIR/raw_counts.h5ad
# target/executable/normalisation_methods/normalise_by_volume/normalise_by_volume ---engine native --input $OUT_DIR/raw_counts.h5ad --input_volume $OUT_DIR/cell_volumes.h5ad --output $OUT_DIR/norm_counts.h5ad
# target/executable/data_processors/process_dataset/process_dataset ---engine native --input_sc resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad --input_sp resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr --output_sc $OUT_DIR/sc_normalised.h5ad --output_sp $OUT_DIR/spatial_processed.zarr
# target/executable/celltype_annotation_methods/ssam/ssam ---engine native --input $OUT_DIR/norm_counts.h5ad --input_transcripts $OUT_DIR/assigned_transcripts.zarr --input_sc $OUT_DIR/sc_normalised.h5ad --output $OUT_DIR/spatial_with_celltypes.h5ad
# target/executable/expr_correction_methods/gene_efficiency_correction/gene_efficiency_correction ---engine native --input $OUT_DIR/spatial_with_celltypes.h5ad --input_sc $OUT_DIR/sc_normalised.h5ad --output $OUT_DIR/spatial_corrected.h5ad