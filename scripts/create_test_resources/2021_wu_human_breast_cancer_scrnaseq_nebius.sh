#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Matching single-cell reference for the atera human breast-cancer spatial test sample.
# Runs on Nebius via tw launch and publishes to scratch (no local run).
#
# The process_wu_human_breast_cancer_sc workflow does NOT subsample internally, so this
# is a two-step chain (as with the yao mouse-brain test SC, which only gets a small
# reference because its loader subsamples):
#   1. run the workflow -> full processed SC atlas on scratch
#   2. chain the openproblems `subsample` module -> ~400-cell test reference on scratch
publish_dir="/scratch/task_ist_preprocessing/resources_test/common"

full_id="2021_wu_human_breast_cancer_scrnaseq_full"
test_id="2021_wu_human_breast_cancer_scrnaseq"

# ---------------------------------------------------------------------------
# Step 1: process the full Wu 2021 breast-cancer scRNA-seq atlas
# ---------------------------------------------------------------------------
cat > /tmp/params_wu_full.yaml << HERE
param_list:
  - id: $full_id
    dataset_name: Wu 2021 Human Breast Cancer scRNAseq
    dataset_organism: Homo sapiens

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

#tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
#  --revision build/main \
#  --pull-latest \
#  --main-script target/nextflow/datasets/workflows/process_wu_human_breast_cancer_sc/main.nf \
#  --workspace 167877437119966 \
#  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
#  --params-file /tmp/params_wu_full.yaml \
#  --config src/base/labels_nebius.config \
#  --labels test_resources,wu_human_breast_cancer_sc \
#  --wait SUCCEEDED

# ---------------------------------------------------------------------------
# Step 2: subsample the processed atlas to a small test reference (~400 cells).
# Requires step 1 to have finished (see --wait SUCCEEDED above); if your tw CLI
# predates --wait, drop it and launch this step manually once step 1 completes.
# ---------------------------------------------------------------------------
cat > /tmp/params_wu_subsample.yaml << HERE
param_list:
  - id: $test_id
    input: $publish_dir/$full_id/dataset.h5ad

n_obs: 400
n_vars: 10000
seed: 0
output: "\$id/dataset.h5ad"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/processors/subsample/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_wu_subsample.yaml \
  --config src/base/labels_nebius.config \
  --labels test_resources,subsample

# The small test reference lands under $publish_dir/$test_id/dataset.h5ad.
# To publish it as a test resource afterwards, sync it from scratch to S3, e.g.:
#   aws s3 sync --profile op \
#     "$publish_dir/$test_id" \
#     "s3://openproblems-data/resources_test/common/$test_id" \
#     --delete --dryrun
