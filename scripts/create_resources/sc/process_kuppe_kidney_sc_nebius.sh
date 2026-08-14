#!/bin/bash

# Process the Kuppe kidney snRNA-seq references (cancer + healthy) through the
# generic SC workflow (log_cp -> hvg -> pca -> knn -> extract_uns_metadata),
# turning the standardized raw counts files into file_common_scrnaseq datasets
# (with a 'normalized' layer, HVG, PCA and kNN) usable as input_sc for the
# corresponding process_datasets_*_nebius.sh combine step.
#
# Prerequisite: run scripts/create_resources/sc/standardize_kuppe_kidney_sc.py
# locally on each Flex-snRNA/{Cancer,Healthy}/annotated.h5ad and upload its output
# to the S3 paths referenced below (see that script's docstring).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# standardized (raw, counts-only) SC references produced by standardize_kuppe_kidney_sc.py
raw_dir="s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney"
publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_kuppe_kidney_sc.yaml << HERE
param_list:

  - id: "kuppe_kidney_cancer_sc/annotated"
    input: "$raw_dir/cancer/dataset.h5ad"

  - id: "kuppe_kidney_healthy_sc/annotated"
    input: "$raw_dir/healthy/dataset.h5ad"

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_scrnaseq/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_kuppe_kidney_sc.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,kuppe_kidney_sc
