#!/bin/bash

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --id mouse_brain_combined \
  --input_sc resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad \
  --input_sp resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr \
  --dataset_id brain_rep1 \
  --dataset_name "Brainz" \
  --dataset_url "https://example.com" \
  --dataset_reference "10.0123456/789456123" \
  --dataset_summary "A summary" \
  --dataset_description "A description" \
  --dataset_organism "mus_musculus" \
  --output_sp '$id/output_sp.zarr' \
  --output_sc '$id/output_sc.h5ad' \
  --output_state '$id/state.yaml' \
  --publish_dir temp_output
