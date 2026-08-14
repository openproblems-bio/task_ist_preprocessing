#!/bin/bash

# Process ONLY the Kuppe kidney MERFISH datasets (combine step) — TWO condition-specific pairings.
# Unlike the LTX/MPII combine scripts (one shared SC reference across samples), each Kuppe sample
# is paired with its OWN condition-matched SC reference, because the cancer and healthy releases use
# disjoint cell-type vocabularies AND different gene panels:
#   cancer  : MERFISH 500-gene Pan-Cancer panel  <-> kuppe_kidney_cancer_sc  (Flex snRNA-seq, Cancer)
#   healthy : MERFISH 140-gene Spapros panel      <-> kuppe_kidney_healthy_sc (Flex snRNA-seq, Healthy)
#
# Reads each spatial input (process_vizgen_merscope loader output) from the local /scratch raw
# folder produced by process_kuppe_kidney_merfish_nebius.sh, and combines it with the matching
# processed SC reference produced by process_kuppe_kidney_sc_nebius.sh.
#
# Prerequisites (both publish to the same /scratch raw folder used below):
#   1. scripts/create_resources/spatial/process_kuppe_kidney_merfish_nebius.sh
#        -> /scratch/.../raw/vizgen_merscope/kuppe_kidney_{cancer,healthy}_merfish/rep1/dataset.zarr
#   2. scripts/create_resources/sc/process_kuppe_kidney_sc_nebius.sh  (log_cp -> hvg -> pca -> knn,
#      i.e. it must carry a 'normalized' layer)
#        -> /scratch/.../raw/kuppe_kidney_{cancer,healthy}_sc/annotated/dataset.h5ad

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# process_kuppe_kidney_merfish_nebius.sh + process_kuppe_kidney_sc_nebius.sh both publish here.
raw_dir='/scratch/task_ist_preprocessing/raw'
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

cat > /tmp/params_kuppe_kidney.yaml << HERE
param_list:

  - id: "kuppe_kidney_cancer_merfish_combined/rep1"
    input_sp: "$raw_dir/vizgen_merscope/kuppe_kidney_cancer_merfish/rep1/dataset.zarr"
    input_sc: "$raw_dir/kuppe_kidney_cancer_sc/annotated/dataset.h5ad"
    dataset_id: "kuppe_kidney_cancer_merfish_combined/rep1"
    dataset_name: "Kidney cancer combined Kuppe MERFISH + Kuppe snRNAseq (Cancer)"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "Kuppe Clear Cell Kidney Carcinoma MERFISH (Vizgen MERSCOPE, 500-gene Pan-Cancer panel) + Kuppe kidney cancer snRNAseq reference"
    dataset_description: "Kuppe lab Clear Cell Kidney Carcinoma FFPE sample (MERFISH, Vizgen MERSCOPE, 500-gene Pan-Cancer panel) paired with the Kuppe kidney cancer Flex snRNA-seq reference. Note: the Pan-Cancer panel does not separate NKT vs T cells or pericytes vs vSMC."
    dataset_organism: "homo_sapiens"

  - id: "kuppe_kidney_healthy_merfish_combined/rep1"
    input_sp: "$raw_dir/vizgen_merscope/kuppe_kidney_healthy_merfish/rep1/dataset.zarr"
    input_sc: "$raw_dir/kuppe_kidney_healthy_sc/annotated/dataset.h5ad"
    dataset_id: "kuppe_kidney_healthy_merfish_combined/rep1"
    dataset_name: "Kidney healthy combined Kuppe MERFISH + Kuppe snRNAseq (Healthy)"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "Kuppe healthy kidney MERFISH (Vizgen MERSCOPE, 140-gene Spapros panel) + Kuppe kidney healthy snRNAseq reference"
    dataset_description: "Kuppe lab healthy kidney FFPE sample (MERFISH, Vizgen MERSCOPE, 140-gene Spapros panel) paired with the Kuppe kidney healthy Flex snRNA-seq reference."
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_kuppe_kidney.yaml "kuppe_kidney"
