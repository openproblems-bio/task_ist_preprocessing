#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/task_ist_preprocessing/results/${RUN_ID}"

cat > /tmp/params_settings.yaml << HERE
default_methods:
  - custom_segmentation
  - basic_transcript_assignment
  - basic_count_aggregation
  - basic_qc_filter
  - alpha_shapes
  - normalize_by_volume
  - ssam
  - no_correction
segmentation_methods:
  - custom_segmentation
  - cellpose
  - binning
  - stardist
  - watershed
transcript_assignment_methods:
  - basic_transcript_assignment
  - baysor
  - clustermap
  - pciseq
  - comseg
  - proseg
count_aggregation_methods:
  - basic_count_aggregation
qc_filtering_methods:
  - basic_qc_filter
volume_calculation_methods:
  - alpha_shapes
normalization_methods:
  - normalize_by_volume
  - normalize_by_counts
  - spanorm
celltype_annotation_methods:
  - ssam
  - tacco
  - moscot
  - mapmycells
  - tangram
  - singler
expression_correction_methods:
  - no_correction
  - gene_efficiency_correction
  - resolvi_correction
method_parameters_yaml: /tmp/method_params.yaml
HERE

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/task_ist_preprocessing/datasets/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# NOTE: this file needs to be made available on the seqera cloud workspace and the 
#       path needs to be added above (method_parameters_yaml)
cat > /tmp/method_params.yaml << HERE
parameters:
  binning:
    default:
      bin_size: 30
    sweep:
      bin_size: [20, 30, 40]
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_ist_preprocessing,full

aws s3 sync \
  s3://openproblems-data/resources/task_ist_preprocessing/results \
  resources/task_ist_preprocessing/results \
  --profile op \
  --dryrun
