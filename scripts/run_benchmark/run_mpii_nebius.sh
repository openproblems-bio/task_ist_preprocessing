#!/bin/bash

# Run the FULL benchmark (all methods) on the MPII + LTX human-lung Xenium datasets.
#
# Same method settings as run_full_nebius.sh; the only difference is that
# `input_states` is globbed to the MPII and LTX combined datasets instead of all
# datasets under the resources dir. These are published by
# scripts/create_resources/combine/process_datasets_mpii_nebius.sh (id
# "mpii_human_lung_xenium_combined/978_reg1") and process_datasets_ltx_nebius.sh
# (ids "ltx_human_lung_xenium_combined/{B_AK-16-12103,L_AK14_14254}"), so their
# state.yaml files live under
# $resources_s3/{mpii,ltx}_human_lung_xenium_combined/**/state.yaml.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_mpii_ltx"

# restrict the benchmark to the MPII and LTX combined datasets only.
# The brace pattern is expanded by Nextflow's glob (not the shell — it stays
# literal inside the heredoc below), matching both the MPII sample and the two
# LTX samples under their respective combined-dataset folders.
dataset_glob="$resources_s3/{mpii_human_lung_xenium_combined,ltx_human_lung_xenium_combined}/**/state.yaml"

cat > /tmp/params_settings_mpii.yaml << HERE
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
  - rctd
expression_correction_methods:
  - no_correction
  - resolvi_correction
  - split
gene_efficiency_correction_methods:
  - no_correction
  - gene_efficiency_correction

#method_parameters_yaml: /tmp/method_params.yaml
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_mpii_benchmark.yaml << HERE
input_states: $dataset_glob
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_mpii.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_mpii_benchmark.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,mpii,ltx
