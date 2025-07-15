#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"


cat > /tmp/params.yaml << HERE
param_list:
  - id: zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc
        
input: "ftp://anonymous@ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad"
keep_files: false 

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_zuani_human_nsclc_sc/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,zuani_human_nsclc_sc

#aws s3 sync \
#  s3://openproblems-data/resources/datasets/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc \
#  resources/datasets/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc
