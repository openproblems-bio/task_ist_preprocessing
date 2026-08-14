#!/bin/bash

# FULL segger run on a few real datasets — validation for the transcript-overlap fix.
#
# PURPOSE: confirm whether the process_dataset crop fix (commit a459afca5, the
# transform-agnostic crop_points_by_global_xy that now crops transcripts to the SAME global
# box as the labels) actually eliminates the out-of-bounds (OOB) transcripts that drove
# segger's `in_bounds` workaround. If it does, segger's OOB-exclusion + row_index remap can be
# simplified to a plain edge clamp (like baysor/proseg). This run is the evidence for that
# decision — DO NOT remove the workaround until this run confirms OOB ~ 0 on real data.
#
# HOW TO READ THE RESULT: segger's script.py prints, per dataset, a line like
#     "<n_oob>/<total> transcripts fall outside the HxW label image; excluding them ..."
# Pull it from each segger task's .command.log (Seqera UI or the work dir). Interpretation:
#   * n_oob ~ 0 (only a handful, edge rounding)  -> crop fix works; safe to simplify the guard.
#   * n_oob still large (e.g. tens of %)         -> either these datasets were NOT re-processed
#                                                   with a459afca5, or the fix is incomplete;
#                                                   KEEP the workaround.
# Also: a completed segger task (no empty-`bd` mid-training crash) is itself part of the signal.
#
# PREREQUISITE (important): the datasets under $resources_s3 must have been RE-PROCESSED with
# the fixed process_dataset (a459afca5). This run consumes already-processed state.yaml files;
# it does NOT re-run process_datasets. If /scratch still holds pre-fix datasets, the OOB count
# reflects the OLD crop, not the fix.
#
# Scope: segger is the ONLY non-default method (the workflow allows at most one non-default
# per pipeline), so per dataset this produces the all-default backbone + the segger variant.
# segger is GPU-only; the `gpuh100` label in its config + src/base/labels_nebius.config
# (runAsUser:0 + /dev/shm volume) pin it to the GPU node group — no GPU-specific change here
# (same mechanism as run_test_segger_nebius.sh / run_segger_fastreseg_comseg_moscot_nebius.sh).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Full datasets live on the Nebius shared /scratch mount (not the small S3 test set).
resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_segger_full_validation"

# --- Datasets to validate on (edit this list to pick 2-3) -----------------------------------
# Chosen to exercise the crop across BOTH transform types the overlap bug hit differently:
#   * 2023_10x_mouse_brain_xenium_combined/rep3  -> Xenium / Scale; the exact dataset segger's
#                                                   NOTES measured 7.6M/19.7M (~38%) OOB on.
#   * 2023_10x_human_breast_cancer_xenium_combined -> a second, larger Xenium (different panel).
#   * 2022_vizgen_human_breast_cancer_merfish_combined/rep1 -> MERFISH / Affine (the OTHER crop
#                                                   regression path; the lung merfish is avoided
#                                                   because it has all-zero labels -> segger's
#                                                   empty-segmentation guard would fail fast).
# IDs are the dataset dir names under $resources_s3 (state.yaml lives at <id>/state.yaml).
datasets=(
  "2023_10x_mouse_brain_xenium_combined/rep3"
  "2023_10x_human_breast_cancer_xenium_combined"
  "2022_vizgen_human_breast_cancer_merfish_combined/rep1"
)

# Build a findStates `--filter` regex (matched against the FULL state-file path with
# .matches(), so anchor with .* and end at .../state.yaml). IDs contain only [A-Za-z0-9_/],
# no regex metacharacters, so they need no escaping; only the literal dot in state.yaml does.
filter_ids=$(IFS='|'; echo "${datasets[*]}")
filter_regex=".*/datasets/(${filter_ids})/state\\.yaml"
echo "Filtering datasets with regex: $filter_regex"

cat > /tmp/params_settings_segger_full.yaml << HERE
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
  - segger
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
# `filter` restricts the broad **/state.yaml glob to just the datasets listed above.
cat > /tmp/params_segger_full.yaml << HERE
input_states: $resources_s3/**/state.yaml
filter: '$filter_regex'
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_segger_full.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_segger_full.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,full,segger,validation
