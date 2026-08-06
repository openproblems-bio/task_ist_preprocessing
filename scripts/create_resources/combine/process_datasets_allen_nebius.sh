#!/bin/bash

# Process ONLY the Allen Brain Cell Atlas (ABCA) MERFISH datasets (combine step) — FOUR
# whole-brain sections, ALL paired with the SAME mouse-brain SC reference (like the LTX/MPII
# combine scripts, and unlike Kuppe which uses condition-specific references):
#   mouse1_coronal, mouse2_coronal, mouse3_sagittal, mouse4_sagittal
#     <-> allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2 (ABCA 2023 Yao 10xv2 atlas)
#
# Reads each spatial input (process_allen_brain_cell_atlas_merfish loader output) and the SC
# reference from the local /scratch raw folder, and writes the combined datasets to the
# /scratch datasets folder (same layout the other *_nebius.sh combine scripts use).
#
# Prerequisites (both publish to the same /scratch raw folder used below):
#   1. scripts/create_resources/spatial/process_allen_brain_cell_atlas_merfish_nebius.sh
#        -> /scratch/.../raw/allen_brain_cell_atlas_merfish/mouse{1,2,3,4}_*/rep1/dataset.zarr
#   2. scripts/create_resources/sc/process_allen_brain_cell_atlas_brain_nebius.sh  (log_cp -> hvg
#      -> pca -> knn, i.e. it must carry a 'normalized' layer)
#        -> /scratch/.../raw/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# process_allen_brain_cell_atlas_merfish_nebius.sh + process_allen_brain_cell_atlas_brain_nebius.sh
# both publish here.
raw_dir='/scratch/task_ist_preprocessing/raw'
publish_dir='/scratch/task_ist_preprocessing/datasets'

# shared SC reference for all four sections (ABCA 2023 Yao whole-brain 10xv2 atlas)
sc_ref="$raw_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"

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

cat > /tmp/params_allen.yaml << HERE
param_list:

  - id: "allen_brain_cell_atlas_merfish_combined/mouse1_coronal/rep1"
    input_sp: "$raw_dir/allen_brain_cell_atlas_merfish/mouse1_coronal/rep1/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "allen_brain_cell_atlas_merfish_combined/mouse1_coronal/rep1"
    dataset_name: "Mouse brain combined ABCA MERFISH mouse1 coronal + 2023 Yao scRNAseq"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_reference: "10.1038/s41586-023-06812-z"
    dataset_summary: "Allen Brain Cell Atlas whole-brain MERFISH mouse 1 (coronal, ~1100-gene panel) + 2023 Yao mouse-brain scRNAseq"
    dataset_description: "Brain-wide MERFISH spatial transcriptomics (Zhuang lab, Allen Brain Cell Atlas), mouse 1 coronal section, paired with the ABCA 2023 Yao whole-mouse-brain 10xv2 scRNAseq reference."
    dataset_organism: "mus_musculus"

  - id: "allen_brain_cell_atlas_merfish_combined/mouse2_coronal/rep1"
    input_sp: "$raw_dir/allen_brain_cell_atlas_merfish/mouse2_coronal/rep1/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "allen_brain_cell_atlas_merfish_combined/mouse2_coronal/rep1"
    dataset_name: "Mouse brain combined ABCA MERFISH mouse2 coronal + 2023 Yao scRNAseq"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_reference: "10.1038/s41586-023-06812-z"
    dataset_summary: "Allen Brain Cell Atlas whole-brain MERFISH mouse 2 (coronal, ~1100-gene panel) + 2023 Yao mouse-brain scRNAseq"
    dataset_description: "Brain-wide MERFISH spatial transcriptomics (Zhuang lab, Allen Brain Cell Atlas), mouse 2 coronal section, paired with the ABCA 2023 Yao whole-mouse-brain 10xv2 scRNAseq reference."
    dataset_organism: "mus_musculus"

  - id: "allen_brain_cell_atlas_merfish_combined/mouse3_sagittal/rep1"
    input_sp: "$raw_dir/allen_brain_cell_atlas_merfish/mouse3_sagittal/rep1/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "allen_brain_cell_atlas_merfish_combined/mouse3_sagittal/rep1"
    dataset_name: "Mouse brain combined ABCA MERFISH mouse3 sagittal + 2023 Yao scRNAseq"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_reference: "10.1038/s41586-023-06812-z"
    dataset_summary: "Allen Brain Cell Atlas whole-brain MERFISH mouse 3 (sagittal, ~1100-gene panel) + 2023 Yao mouse-brain scRNAseq"
    dataset_description: "Brain-wide MERFISH spatial transcriptomics (Zhuang lab, Allen Brain Cell Atlas), mouse 3 sagittal section, paired with the ABCA 2023 Yao whole-mouse-brain 10xv2 scRNAseq reference."
    dataset_organism: "mus_musculus"

  - id: "allen_brain_cell_atlas_merfish_combined/mouse4_sagittal/rep1"
    input_sp: "$raw_dir/allen_brain_cell_atlas_merfish/mouse4_sagittal/rep1/dataset.zarr"
    input_sc: "$sc_ref"
    dataset_id: "allen_brain_cell_atlas_merfish_combined/mouse4_sagittal/rep1"
    dataset_name: "Mouse brain combined ABCA MERFISH mouse4 sagittal + 2023 Yao scRNAseq"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_reference: "10.1038/s41586-023-06812-z"
    dataset_summary: "Allen Brain Cell Atlas whole-brain MERFISH mouse 4 (sagittal, ~1100-gene panel) + 2023 Yao mouse-brain scRNAseq"
    dataset_description: "Brain-wide MERFISH spatial transcriptomics (Zhuang lab, Allen Brain Cell Atlas), mouse 4 sagittal section, paired with the ABCA 2023 Yao whole-mouse-brain 10xv2 scRNAseq reference."
    dataset_organism: "mus_musculus"

# ABCA whole-brain images are enormous (~83k x 102k px) and the tissue sits off-centre,
# so an image-centred crop can miss it (mouse1: 0 of 42M transcripts). Centre the crop
# on the transcript density instead (opt-in flag; default is image-centred).
tissue_centered_crop: true

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_allen.yaml "allen_brain_cell_atlas_merfish"
