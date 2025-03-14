#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Create local output directory
output_dir="resources/datasets"

if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

cat > /tmp/params.yaml << HERE
param_list:
  - id: wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc
    cancer_subtypes:
      - HER2+
      - TNBC
      - ER+

keep_files: false 

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$output_dir"
HERE

# Run nextflow workflow locally
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_wu_human_breast_cancer_sc/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  -c common/nextflow_helpers/labels_ci.config