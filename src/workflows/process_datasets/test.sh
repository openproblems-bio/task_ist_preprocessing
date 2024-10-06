#!/bin/bash

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --id mouse_brain_combined \
  --input_sc resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad \
  --input_sp resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr \
  --output_sp '$id/output_sp.zarr' \
  --output_sc '$id/output_sc.h5ad' \
  --publish_dir temp_output

# created files:
#   output/my_dataset_id0.process_datasets.output_sp.zarr
#   output/my_dataset_id0.process_datasets.output_sc.h5ad
#   output/my_dataset_id1.process_datasets.output_sp.zarr
#   output/my_dataset_id1.process_datasets.output_sc.h5ad

# created files:
#   output/my_dataset_id0/output_sp.zarr
#   output/my_dataset_id0/output_sc.h5ad
#   output/my_dataset_id1/output_sp.zarr
#   output/my_dataset_id1/output_sc.h5ad
