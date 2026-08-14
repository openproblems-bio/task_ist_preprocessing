#!/bin/bash

# Run the FULL benchmark (all methods) on the Allen Brain Cell Atlas (ABCA) MERFISH
# datasets ONLY.
#
# Same method settings as run_full_nebius.sh; the only difference is that
# `input_states` is globbed to the ABCA MERFISH combined datasets instead of all
# datasets under the resources dir. These are published by
# scripts/create_resources/combine/process_datasets_allen_nebius.sh (ids
# "allen_brain_cell_atlas_merfish_combined/mouse{1,2,3,4}_{coronal,sagittal}/rep1"),
# so their state.yaml files live under
# $resources_s3/allen_brain_cell_atlas_merfish_combined/**/state.yaml.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir_s3="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_allen_merfish"

# restrict the benchmark to the ABCA MERFISH combined datasets only.
# The ** is expanded by Nextflow's glob (not the shell — it stays literal inside
# the heredoc below), matching all four mouse sections (mouse1_coronal,
# mouse2_coronal, mouse3_sagittal, mouse4_sagittal) under the combined folder.
dataset_glob="$resources_s3/allen_brain_cell_atlas_merfish_combined/**/state.yaml"

cat > /tmp/params_settings_allen.yaml << HERE
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
  - no_correction
  - resolvi_correction
  - split
gene_efficiency_correction_methods:
  - no_correction
#  - gene_efficiency_correction

#method_parameters_yaml: /tmp/method_params.yaml
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_allen_benchmark.yaml << HERE
input_states: $dataset_glob
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_allen.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir_s3"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_allen_benchmark.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,allen_merfish
