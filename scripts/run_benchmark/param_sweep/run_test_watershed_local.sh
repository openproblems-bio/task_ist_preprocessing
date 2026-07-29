#!/bin/bash

# Test run: all default methods + watershed segmentation, with a parameter sweep
# over the Tier-1 optimization levers identified for classic marker-controlled
# watershed. See src/methods_segmentation/watershed/NOTES.md ("Optimization /
# tuning").
#
# watershed is a CPU method (no GPU needed): the txsim classic pipeline
# normalize -> contrast -> blur -> threshold -> post-process -> distance transform
# -> local-maxima markers -> watershed -> background filter.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

echo "Running benchmark on test data (defaults + watershed sweep)"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="testrun_watershed_$(date +%Y-%m-%d_%H-%M-%S)"
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
  - watershed
  # - cellpose
  # - cellposev4
  # - binning
  # - stardist
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
method_parameters_yaml: $REPO_ROOT/scripts/run_benchmark/watershed_params.yaml
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

# The watershed parameter sweep lives in a committed file (single source of truth,
# shared with run_test_watershed_nebius.sh): scripts/run_benchmark/watershed_params.yaml
# It defines a `default:` variant + one-arg-at-a-time `sweep:` variants (a "star"
# around the default, not a grid) => 12 watershed variants. Edit that file to
# change the sweep. Referenced via $REPO_ROOT above so local Nextflow reads it directly.

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file /tmp/params.yaml
