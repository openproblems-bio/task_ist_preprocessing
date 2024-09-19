#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

input_dir="s3://openproblems-data/resources/datasets"
publish_dir="s3://openproblems-data/resources/task_ist_preprocessing/datasets"


cat > /tmp/params.yaml << HERE
param_list:

  - id: "mouse_brain_combined/rep1"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep1/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch openproblems-bio/task_ist_preprocessing \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,10x_xenium
