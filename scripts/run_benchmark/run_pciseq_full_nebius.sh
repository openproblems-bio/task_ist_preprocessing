#!/bin/bash

# FULL pciSeq run across ALL datasets on Nebius — pciSeq isolated as the single
# non-default method.
#
# PURPOSE: run pciSeq transcript assignment on every processed dataset on its own, rather
# than as one variant inside the big star-matrix of run_full_nebius.sh. Isolating it makes
# it easy to (re)launch / watch / triage just pciSeq (e.g. after the empty-segmentation guard
# + re-processing fix) without re-running the whole matrix.
#
# Scope: pciSeq is the ONLY non-default method (the workflow allows at most one non-default
# per pipeline), so per dataset this produces the all-default backbone + the pciSeq variant.
# Runs on ALL datasets under $resources_s3 (the broad **/state.yaml glob, no `filter`).
#
# PREREQUISITE: the datasets under $resources_s3 must already be processed (this run consumes
# existing state.yaml files; it does NOT re-run process_datasets).
#
# COST NOTE: pciSeq is CPU (its config `label: [veryhightime, midcpu, highmem]`) — no GPU
# nodes. One pciSeq task per dataset, each up to the veryhightime (24h) ceiling on large panels.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Full datasets live on the Nebius shared /scratch mount (not the small S3 test set).
resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_pciseq_full"

cat > /tmp/params_settings_pciseq_full.yaml << HERE
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
transcript_assignment_methods:
  - basic_transcript_assignment
  - pciseq
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
expression_correction_methods:
  - no_correction
gene_efficiency_correction_methods:
  - no_correction
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this).
cat > /tmp/params_pciseq_full.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_pciseq_full.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_pciseq_full.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,full,pciseq
