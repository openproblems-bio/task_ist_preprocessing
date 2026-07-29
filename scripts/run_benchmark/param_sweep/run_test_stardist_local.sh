#!/bin/bash

# Test run: all default methods + StarDist2D (stardist) segmentation, with a
# parameter sweep over the optimization levers identified for StarDist.
# See src/methods_segmentation/stardist/NOTES.md ("Optimization / tuning").
#
# NOTE: stardist runs on TensorFlow via a GPU-capable base image
# (openproblems/base_tensorflow_nvidia). It falls back to CPU if no GPU is present
# (slower, but far lighter than cellposev4's SAM ViT). Run on a CUDA-capable Docker
# host for the full sweep.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

echo "Running benchmark on test data (defaults + stardist sweep)"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="testrun_stardist_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="temp/results/${RUN_ID}"

cat > /tmp/params_settings.yaml << HERE
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
  - stardist
  # - cellpose
  # - cellposev4
  # - binning
  # - watershed
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
  - ssam
expression_correction_methods:
  - no_correction
gene_efficiency_correction_methods:
  - no_correction
method_parameters_yaml: $REPO_ROOT/scripts/run_benchmark/stardist_params.yaml
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry auto` for this)
cat > /tmp/params.yaml << HERE
input_states: resources_test/task_ist_preprocessing/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# The stardist parameter sweep lives in a committed file (single source of truth,
# shared with run_test_stardist_nebius.sh): scripts/run_benchmark/stardist_params.yaml
# It defines a `default:` variant + one-arg-at-a-time `sweep:` variants (a "star"
# around the default, not a grid) => 13 stardist variants. Edit that file to change
# the sweep. Referenced via $REPO_ROOT above so local Nextflow reads it directly.

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file /tmp/params.yaml
