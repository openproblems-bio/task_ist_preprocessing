#!/bin/bash

# TODO: The param_list metadata was mostly infered with chatGPT from the create resources scripts.
#       Double check if everything's correct.


# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

input_dir="s3://openproblems-data/resources/datasets"
#publish_dir="s3://openproblems-data/resources/task_ist_preprocessing/datasets"
publish_dir='/scratch/task_ist_preprocessing/datasets'

cat > /tmp/params_bruker.yaml << HERE
param_list:

  - id: "bruker_mouse_brain_cosmx_combined/rep1"
    input_sp: "$input_dir/bruker_cosmx/bruker_mouse_brain_cosmx/rep1/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_id: "bruker_mouse_brain_cosmx_combined/rep1"
    dataset_name: "Mouse brain combined Bruker CosMx rep1 2023 Yao scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-smi-mouse-brain-ffpe-dataset/"
    dataset_reference: "10.1038/s41586-023-06812-z"
    dataset_summary: "Bruker CosMx Mouse Brain + ABCA Mouse Brain scRNAseq"
    dataset_description: "Bruker CosMx Mouse Brain + ABCA Mouse Brain scRNAseq"
    dataset_organism: "mus_musculus"

  - id: "bruker_human_liver_cosmx_combined"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_liver_cosmx/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Andrews_human_liver_sc/dataset.h5ad"
    dataset_id: "bruker_human_liver_cosmx_combined"
    dataset_name: "Human liver combined Bruker CosMx 2022 Andrews scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/"
    dataset_reference: "https://doi.org/10.1002/hep4.1854"
    dataset_summary: "Bruker CosMx Human Liver + 2022 Andrews scRNAseq"
    dataset_description: "Bruker CosMx Human Liver + 2022 Andrews scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_liver_cancer_cosmx_combined"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_liver_cancer_cosmx/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Lu_human_liver_cancer_sc/dataset.h5ad"
    dataset_id: "bruker_human_liver_cancer_cosmx_combined"
    dataset_name: "Human liver cancer combined Bruker CosMx 2022 Lu scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/"
    dataset_reference: "https://doi.org/10.1038/s41467-022-32283-3"
    dataset_summary: "Bruker CosMx Human Liver Cancer + 2022 Lu scRNAseq"
    dataset_description: "Bruker CosMx Human Liver Cancer + 2022 Lu scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep1"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep1/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep1"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung5 rep1 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep1 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep1 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep2"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep2/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep2"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung5 rep2 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep2 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep2 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep3"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep3/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung5_rep3"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung5 rep3 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep3 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep3 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung6"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung6/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung6"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung6 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung6 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung6 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung9_rep1"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung9_rep1/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung9_rep1"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung9 rep1 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung9 Rep1 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung9 Rep1 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung9_rep2"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung9_rep2/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung9_rep2"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung9 rep2 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung9 Rep2 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung9 Rep2 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung12"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung12/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung12"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung12 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung12 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung12 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "bruker_human_lung_cancer_cosmx_combined/lung13"
    input_sp: "$input_dir/bruker_cosmx/bruker_human_lung_cancer_cosmx/lung13/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "bruker_human_lung_cancer_cosmx_combined/lung13"
    dataset_name: "Human lung cancer combined Bruker CosMx Lung13 2024 Zuani scRNAseq"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung13 + 2024 Zuani scRNAseq"
    dataset_description: "Bruker CosMx Human Lung Cancer Lung13 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_bruker.yaml \
  --config src/base/labels_nebius.config \
  --labels "task_ist_preprocessing,process_datasets,bruker"