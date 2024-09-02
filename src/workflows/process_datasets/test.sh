#!/bin/bash


cat > /tmp/params.yaml <<EOF
param_list:
  - id: my_dataset_id0
    input_sp: path/to/my_dataset0/dataset.h5ad
    input_sc: path/to/my_dataset0/dataset.csv
  - id: my_dataset_id1
    input_sp: path/to/my_dataset1/dataset.h5ad
    input_sc: path/to/my_dataset1/dataset.csv
reference: my_reference.csv
EOF

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -params-file /tmp/params.yaml \
  --publish_dir output \
  --output_sp '$id/output_sp.zarr' \
  --output_sc '$id/output_sc.h5ad'

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
