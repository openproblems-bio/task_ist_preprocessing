#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/task_ist_preprocessing/results/${RUN_ID}"

# input_dir="s3://openproblems-data/resources/task_ist_preprocessing/datasets"
# cat > /tmp/params.yaml << HERE
# param_list:

#   - id: "mouse_brain_combined/rep1"
#     input_sp: "$input_dir/mouse_brain_combined/rep1/output_sp.zarr"
#     input_sc: "$input_dir/mouse_brain_combined/rep1/output_sc.h5ad"

# output_sc: "\$id/output_sc.h5ad"
# output_sp: "\$id/output_sp.zarr"
# output_state: "\$id/state.yaml"
# publish_dir: "$publish_dir"
# HERE

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/task_ist_preprocessing/datasets/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_ist_preprocessing,full

aws s3 sync \
  s3://openproblems-data/resources/task_ist_preprocessing/results \
  resources/task_ist_preprocessing/results \
  --profile op \
  --dryrun
