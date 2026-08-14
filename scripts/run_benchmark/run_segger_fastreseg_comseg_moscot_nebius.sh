#!/bin/bash

# Full (ALL datasets) run for the segger, fastreseg, comseg (transcript-assignment)
# and moscot (cell-type-annotation) methods. Runs over the full datasets on the
# Nebius shared /scratch mount, using the DEFAULT method at every stage, and adds:
#   * segger, fastreseg, comseg to the transcript-assignment stage
#   * moscot to the cell-type-annotation stage
#
# The benchmark workflow enforces "at most one non-default method per pipeline"
# (see expandChannelWithParameterSets in src/workflows/run_benchmark/main.nf),
# i.e. it varies exactly one stage at a time from the defaults. So, per dataset,
# this launch produces these runs:
#   - all-default backbone
#   - + segger      (transcript assignment; rest default, incl. tacco)
#   - + fastreseg   (transcript assignment; rest default, incl. tacco)
#   - + comseg      (transcript assignment; rest default, incl. tacco) -> 
#   - + moscot      (cell-type annotation; rest default, incl. basic_transcript_assignment)
# There is no segger+moscot etc. combination — that would be two non-default steps.
#
# Resource routing (all handled by src/base/labels_nebius.config, one launch):
#   * segger    -> GPU  (L40S node group; [hightime, midcpu, highmem, gpu])
#   * moscot    -> GPU  (H100 node group; [hightime, midcpu, veryhighmem, gpuh100])
#   * fastreseg -> CPU  ([midtime, midcpu, midmem])
#   * comseg    -> CPU  ([hightime, midcpu, highmem])

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Full datasets live on the Nebius shared /scratch mount (not the small S3 test set).
#resources_s3=s3://openproblems-data/resources_test/task_ist_preprocessing
resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_segger_fastreseg_comseg_moscot"

cat > /tmp/params_settings_segger_fastreseg_comseg_moscot.yaml << HERE
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
#  - moscot
expression_correction_methods:
  - no_correction
gene_efficiency_correction_methods:
  - no_correction
  # - gene_efficiency_correction
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_segger_fastreseg_comseg_moscot.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_segger_fastreseg_comseg_moscot.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_segger_fastreseg_comseg_moscot.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,full,segger,fastreseg,comseg,moscot
