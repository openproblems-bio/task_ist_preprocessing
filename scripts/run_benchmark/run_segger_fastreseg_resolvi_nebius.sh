#!/bin/bash

# segger + fastreseg transcript-assignment run across ALL datasets on Nebius, with resolVI as the
# DEFAULT expression-correction method (the sibling run_segger_fastreseg_comseg_moscot_nebius.sh
# uses `no_correction`).
#
# Same shape as that script: the benchmark enforces "at most one non-default method per
# pipeline" (expandChannelWithParameterSets in src/workflows/run_benchmark/main.nf), so per
# dataset this launch produces:
#   - all-default backbone   (basic_transcript_assignment; rest default) + resolVI
#   - + segger               (transcript assignment; rest default)       + resolVI
#   - + fastreseg            (transcript assignment; rest default)       + resolVI
# resolVI is applied to EVERY one of these because it is the always-run backbone correction.
#
# HOW THE DEFAULT IS SET (same idiom as run_full_resolvi_nebius.sh / run_rctd_split_nebius.sh):
#   `resolvi_correction` is in `default_methods` AND is the sole entry of
#   `expression_correction_methods`, so it is the always-run backbone correction. The expr-corr
#   `no_correction` component is not selected at that stage, so it never runs there;
#   `no_correction` stays in `default_methods` only as the gene-efficiency default (a different
#   component that shares the config name).
#
# Resource routing (all via src/base/labels_nebius.config, one launch):
#   * segger    -> GPU  (L40S node group; [hightime, midcpu, highmem, gpu])
#   * fastreseg -> CPU  ([midtime, midcpu, midmem])
#   * resolvi   -> GPU  ([hightime, highcpu, highmem, gpu]) -- runs on EVERY branch (default),
#                       so the segger branch is GPU at BOTH transcript-assignment and correction.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Full datasets live on the Nebius shared /scratch mount (not the small S3 test set).
resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_segger_fastreseg_resolvi"

cat > /tmp/params_settings_segger_fastreseg_resolvi.yaml << HERE
default_methods:
  - custom_segmentation
  - basic_transcript_assignment
  - basic_count_aggregation
  - basic_qc_filter
  - alpha_shapes
  - normalize_by_volume
  - tacco
  - resolvi_correction
  - no_correction
segmentation_methods:
  - custom_segmentation
transcript_assignment_methods:
  - basic_transcript_assignment
  - segger
  - fastreseg
# - comseg
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
  - resolvi_correction
gene_efficiency_correction_methods:
  - no_correction
  # - gene_efficiency_correction
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_segger_fastreseg_resolvi.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_segger_fastreseg_resolvi.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_segger_fastreseg_resolvi.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,segger,fastreseg,resolvi
