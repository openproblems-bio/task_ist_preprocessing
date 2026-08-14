#!/bin/bash

# Process ONLY the custom LTX human-lung Xenium dataset (combine step) — TWO samples.
# Reads each spatial input (process_tenx_xenium loader output) from the local
# /scratch raw folder produced by process_ltx_xenium_nebius.sh, and combines each
# with the (single, shared) LTX human lung scRNAseq reference.
#
# The scRNAseq reference must already be standardized (scripts/create_resources/
# sc/standardize_ltx_sc.py) AND processed (log_cp -> hvg -> pca -> knn, i.e. it
# must carry a 'normalized' layer) before it is used here — run
# scripts/create_resources/sc/process_ltx_sc_nebius.sh first.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Xenium loader output (process_ltx_xenium_nebius.sh publishes here).
raw_dir='/scratch/task_ist_preprocessing/raw'
# scRNAseq reference is the process_ltx_sc_nebius.sh output (log_cp/hvg/pca/knn).
sc_ref="/scratch/task_ist_preprocessing/raw/ltx_human_lung_sc/ltx_full_annotation_2026_04_25/dataset.h5ad"
publish_dir='/scratch/task_ist_preprocessing/datasets'

launch_batch() {
  local params_file="$1"
  local label="$2"
  tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/process_datasets/main.nf \
    --workspace 167877437119966 \
    --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
    --params-file "$params_file" \
    --config src/base/labels_nebius.config \
    --labels "task_ist_preprocessing,process_datasets,$label"
}

cat > /tmp/params_ltx.yaml << HERE
param_list:

  - id: "ltx_human_lung_xenium_combined/B_AK-16-12103"
    input_sp: "$raw_dir/ltx_human_lung_xenium/B_AK-16-12103/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "ltx_human_lung_xenium_combined/B_AK-16-12103"
    dataset_name: "Human lung combined LTX Xenium B_AK-16-12103 + LTX scRNAseq"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "LTX Xenium human lung (sample B_AK-16-12103) + LTX human lung scRNAseq reference"
    dataset_description: "LTX Xenium human lung (sample B_AK-16-12103) + LTX human lung scRNAseq reference (ltx_full_annotation_2026_04_25)"
    dataset_organism: "homo_sapiens"

  - id: "ltx_human_lung_xenium_combined/L_AK14_14254"
    input_sp: "$raw_dir/ltx_human_lung_xenium/L_AK14_14254/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "ltx_human_lung_xenium_combined/L_AK14_14254"
    dataset_name: "Human lung combined LTX Xenium L_AK14_14254 + LTX scRNAseq"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "LTX Xenium human lung (sample L_AK14_14254) + LTX human lung scRNAseq reference"
    dataset_description: "LTX Xenium human lung (sample L_AK14_14254) + LTX human lung scRNAseq reference (ltx_full_annotation_2026_04_25)"
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_ltx.yaml "ltx"
