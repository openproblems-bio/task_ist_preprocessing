#!/bin/bash

# Process ONLY the MPII custom human-lung Xenium dataset (combine step).
# Reads the spatial input (process_tenx_xenium loader output) from the local
# /scratch raw folder produced by process_mpii_xenium_nebius.sh, and combines it
# with the MPII human lung scRNAseq reference.
#
# The scRNAseq reference must already be standardized (scripts/create_resources/
# sc/standardize_mpii_human_lung_sc.py) AND processed (log_cp -> hvg -> pca ->
# knn, i.e. it must carry a 'normalized' layer) before it is used here.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Xenium loader output (process_mpii_xenium_nebius.sh publishes here).
raw_dir='/scratch/task_ist_preprocessing/raw'
# scRNAseq reference is not a loader output; it lives in the S3 datasets bucket.
# Update this to wherever the standardized+processed reference was uploaded.
sc_ref="/scratch/task_ist_preprocessing/raw/mpii_human_lung_sc/current_annotation/dataset.h5ad"
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

cat > /tmp/params_mpii.yaml << HERE
param_list:

  - id: "mpii_human_lung_xenium_combined/978_reg1"
    input_sp: "$raw_dir/mpii_human_lung_xenium/978_reg1/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "mpii_human_lung_xenium_combined/978_reg1"
    dataset_name: "Human lung combined MPII Xenium 978 reg1 + MPII scRNAseq"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "MPII Xenium human lung (sample 978, region 1) + MPII human lung scRNAseq reference"
    dataset_description: "MPII Xenium human lung (sample 978, region 1) + MPII human lung scRNAseq reference"
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_mpii.yaml "mpii"
