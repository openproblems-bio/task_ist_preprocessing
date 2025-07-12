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
  - id: process_zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc
        
input: "ftp://anonymous@ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad"
keep_files: false 

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$output_dir"
HERE

# Run nextflow workflow locally
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_zuani_human_nsclc_sc/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  -c common/nextflow_helpers/labels_ci.config