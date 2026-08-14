#!/bin/bash

# Process the LTX human-lung scRNAseq reference through the generic SC workflow
# (log_cp -> hvg -> pca -> knn -> extract_uns_metadata), turning the standardized
# raw counts file into a file_common_scrnaseq dataset (with a 'normalized' layer,
# HVG, PCA and kNN) usable as input_sc for process_datasets_ltx_nebius.sh.
#
# Prerequisite: run scripts/create_resources/sc/standardize_ltx_sc.py locally and
# upload its output to the S3 path referenced by `input` below (the same path the
# standardize script's docstring uploads to).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# standardized (raw, counts-only) SC reference produced by standardize_ltx_sc.py
input_h5ad="s3://openproblems-data/resources/raw_data/txSim_custom/ltx/dataset.h5ad"
publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_ltx_sc.yaml << HERE
param_list:

  - id: "ltx_human_lung_sc/ltx_full_annotation_2026_04_25"
    input: "$input_h5ad"

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
  --params-file /tmp/params_ltx_sc.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,ltx_sc
