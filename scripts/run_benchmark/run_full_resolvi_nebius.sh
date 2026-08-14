#!/bin/bash

# Run the FULL benchmark matrix across ALL datasets on Nebius, but with resolVI as the
# DEFAULT expression-correction method (run_full_nebius.sh uses `no_correction`).
#
# Modeled on run_full_nebius.sh: every per-stage sweep list is IDENTICAL, so the same
# star-matrix of segmentation / transcript-assignment / normalization / cell-type-annotation
# variants is explored. The ONLY change is the expression-correction step:
#   run_full_nebius.sh : expr-corr default = no_correction     (correction OFF by default)
#   this script        : expr-corr default = resolvi_correction (resolVI applied to every pipeline)
#
# HOW THE DEFAULT IS SET (same idiom as run_rctd_split_nebius.sh):
#   `resolvi_correction` is placed in `default_methods` AND is the sole entry of
#   `expression_correction_methods`, so it is the always-run backbone correction that every
#   swept pipeline chains through. The expr-corr `no_correction` component is simply not
#   selected at that stage, so it never runs there. `no_correction` stays in `default_methods`
#   only as the gene-efficiency default -- a DIFFERENT component that happens to share the
#   config name `no_correction`; it runs only at the gene-efficiency stage, the sole stage
#   that still selects it.
#
# COST NOTE: resolVI is GPU + hightime + highmem (see its config `label: [hightime, highcpu,
# highmem, gpu]`; GPU tasks are routed by src/base/labels_nebius.config, same as run_gpu_nebius.sh).
# As the DEFAULT it now runs on EVERY pipeline in the matrix, so this schedules a LARGE number
# of GPU jobs -- considerably more expensive than run_full or run_gpu (which only sweep resolVI
# against the all-default backbone). Trim the upstream sweep lists if that is too heavy.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_resolvi"

cat > /tmp/params_settings_resolvi.yaml << HERE
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
  - cellpose
  - cellposev4
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
  - segger
  - fastreseg
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
  - rctd
expression_correction_methods:
  - resolvi_correction
gene_efficiency_correction_methods:
  - no_correction
#  - gene_efficiency_correction
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_resolvi.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_resolvi.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_resolvi.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,resolvi
