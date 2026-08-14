#!/bin/bash

# Run a SINGLE-pipeline benchmark across ALL datasets on Nebius:
#   - cell type annotation : RCTD
#   - expression correction: SPLIT
#   - every other step uses its DEFAULT method (no method sweeps)
#
# Modeled on run_full_nebius.sh, but each per-stage method list is a singleton
# so the star-matrix collapses to exactly one pipeline per dataset. For that
# single chain to survive every stage, ALL of its methods must appear in
# `default_methods` (the star-matrix truncates once any non-default method is
# used upstream). `no_correction` is the gene-efficiency default; it only runs
# at the gene-eff stage because that is the only stage that selects it, so it is
# harmless alongside `split` in default_methods.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_rctd_split"

cat > /tmp/params_settings_rctd_split.yaml << HERE
default_methods:
  - custom_segmentation
  - basic_transcript_assignment
  - basic_count_aggregation
  - basic_qc_filter
  - alpha_shapes
  - normalize_by_volume
  - rctd
  - split
  - no_correction
segmentation_methods:
  - custom_segmentation
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
  - rctd
expression_correction_methods:
  - split
gene_efficiency_correction_methods:
  - no_correction
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_rctd_split.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_rctd_split.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_rctd_split.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,rctd,split
