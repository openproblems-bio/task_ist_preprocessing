#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_dir="/scratch/task_ist_preprocessing/datasets"
publish_dir="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_gpu"

cat > /tmp/params_settings_gpu.yaml << HERE
default_methods:
  - custom_segmentation
  - basic_transcript_assignment
  - basic_count_aggregation
  - basic_qc_filter
  - alpha_shapes
  - normalize_by_volume
  - tacco
  - no_correction
segmentation_methods:
  - custom_segmentation
  - cellpose
  - cellposev4
  - stardist
transcript_assignment_methods:
  - basic_transcript_assignment
count_aggregation_methods:
  - basic_count_aggregation
qc_filtering_methods:
  - basic_qc_filter
volume_calculation_methods:
  - alpha_shapes
normalization_methods:
  - normalize_by_volume
celltype_annotation_methods:
  - tacco
  - tangram
  - moscot
expression_correction_methods:
  - no_correction
  - resolvi_correction
gene_efficiency_correction_methods:
  - no_correction
  # - gene_efficiency_correction
HERE

cat > /tmp/params_gpu.yaml << HERE
input_states: $resources_dir/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_gpu.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_gpu.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,gpu