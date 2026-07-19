#!/bin/bash

# Process ONLY the 10x Atera breast cancer dataset.
# Reads the spatial input (Atera loader output) from the local /scratch raw
# folder produced by process_10x_atera_nebius.sh, and combines it with the
# 2021 Wu scRNAseq reference.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Atera loader output (process_10x_atera_nebius.sh publishes here).
raw_dir='/scratch/task_ist_preprocessing/raw'
# scRNAseq reference is not a loader output; it lives in the S3 datasets bucket.
sc_dir="s3://openproblems-data/resources/datasets"
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

cat > /tmp/params_atera.yaml << HERE
param_list:

  - id: "2026_10x_human_breast_cancer_atera_combined"
    input_sp: "$raw_dir/10x_atera/2026_10x_human_breast_cancer_atera/dataset.zarr"
    input_sc: "$sc_dir/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc/dataset.h5ad"
    dataset_id: "2026_10x_human_breast_cancer_atera_combined"
    dataset_name: "Human breast cancer combined 2026 10x Atera WTA 2021 Wu scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer"
    dataset_reference: "https://doi.org/10.1038/s41588-021-00911-1"
    dataset_summary: "Atera WTA FFPE Human Breast Cancer (DCIS Grade 3) + 2021 Wu scRNAseq"
    dataset_description: "Atera WTA FFPE Human Breast Cancer (DCIS Grade 3) + 2021 Wu scRNAseq"
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_atera.yaml "atera"
