#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Local (nextflow run) variant of 2021_wu_human_breast_cancer_scrnaseq_nebius.sh.
# Matching single-cell reference for the atera human breast-cancer spatial test sample.
# Publishes to the local resources_test/ instead of scratch.
#
# The process_wu_human_breast_cancer_sc workflow does NOT subsample internally, so this is
# a two-step chain (as in the nebius script):
#   1. run the workflow -> full processed SC atlas under resources_test/common
#   2. chain the openproblems `subsample` module -> ~400-cell test reference
# `nextflow run` is synchronous, so step 2 runs only after step 1 finishes (no --wait needed).
publish_dir="resources_test/common"

full_id="2021_wu_human_breast_cancer_scrnaseq_full"
test_id="2021_wu_human_breast_cancer_scrnaseq"

# Cap per-process CPU requests to the 8 cores available on this machine. Component
# labels otherwise request up to highcpu=30, which exceeds the machine and the
# Nextflow local executor would then fail to schedule those processes.
cat > /tmp/labels_8cpu.config << 'HERE'
process {
  withLabel: lowcpu  { cpus = 2 }
  withLabel: midcpu  { cpus = 4 }
  withLabel: highcpu { cpus = 8 }
}
executor {
  cpus = 8
}
HERE

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

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_wu_human_breast_cancer_sc/main.nf \
  -profile docker \
  -resume \
  -c /tmp/labels_8cpu.config \
  -params-file /tmp/params_wu_full.yaml

# ---------------------------------------------------------------------------
# Step 2: subsample the processed atlas to a small test reference (~400 cells).
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

nextflow run . \
  -main-script target/nextflow/datasets/processors/subsample/main.nf \
  -profile docker \
  -resume \
  -c /tmp/labels_8cpu.config \
  -params-file /tmp/params_wu_subsample.yaml

# sync the small test reference to s3
aws s3 sync --profile op \
  "resources_test/common/$test_id" \
  "s3://openproblems-data/resources_test/common/$test_id" \
  --delete --dryrun
