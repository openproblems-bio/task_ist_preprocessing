#!/bin/bash

# Process the Ganier human-skin scRNAseq reference on Nebius.
#
# Unlike the custom ltx/mpii SC scripts (which standardize a raw counts file
# locally and push it through the generic process_scrnaseq workflow), Ganier has
# a self-contained dedicated workflow: it downloads via its loader and runs
# hvg -> pca -> knn -> extract_uns_metadata, emitting a full file_common_scrnaseq
# dataset (with a 'normalized' layer, HVG, PCA and kNN) directly.
#
# The dataset source and all metadata defaults live in the workflow config
# (src/datasets/workflows/process_ganier_human_skin_sc/config.vsh.yaml), so only
# the id needs to be specified here.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_ganier_sc.yaml << HERE
param_list:

  - id: scrnaseq_for_ist/2024Ganier_human_skin_sc

keep_files: false

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_ganier_human_skin_sc/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_ganier_sc.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,ganier_human_skin_sc
