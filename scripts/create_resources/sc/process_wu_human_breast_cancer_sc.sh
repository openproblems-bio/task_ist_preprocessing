#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

# Note that the current download script and processing workflow have a specific default parameter set for the given dataset.
# No additional datasets are supported by that component/workflow. Therefore the default parameters are used and don't need 
# to be specified here.

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
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_wu_human_breast_cancer_sc/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,wu_human_breast_cancer_sc

#aws s3 sync \
#  s3://openproblems-data/resources/datasets/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc \
#  resources/datasets/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc
