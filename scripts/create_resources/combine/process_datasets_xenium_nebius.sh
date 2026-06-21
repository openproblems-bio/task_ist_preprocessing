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

cat > /tmp/params_xenium.yaml << HERE
param_list:

  - id: "2023_10x_mouse_brain_xenium_combined/rep1"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep1/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_id: "2023_10x_mouse_brain_xenium_combined/rep1"
    dataset_name: "Mouse brain combined 2023 10x Xenium rep1 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_summary: "Xenium V1 Fresh Frozen Mouse Brain rep1 + ABCA Mouse Brain scRNAseq"
    dataset_description: "Xenium V1 Fresh Frozen Mouse Brain rep1 + ABCA Mouse Brain scRNAseq"
    dataset_organism: "mus_musculus"

  - id: "2023_10x_mouse_brain_xenium_combined/rep2"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep2/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_id: "2023_10x_mouse_brain_xenium_combined/rep2"
    dataset_name: "Mouse brain combined 2023 10x Xenium rep2 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_summary: "Xenium V1 Fresh Frozen Mouse Brain rep2 + ABCA Mouse Brain scRNAseq"
    dataset_description: "Xenium V1 Fresh Frozen Mouse Brain rep2 + ABCA Mouse Brain scRNAseq"
    dataset_organism: "mus_musculus"

  - id: "2023_10x_mouse_brain_xenium_combined/rep3"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep3/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_id: "2023_10x_mouse_brain_xenium_combined/rep3"
    dataset_name: "Mouse brain combined 2023 10x Xenium rep3 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard"
    dataset_summary: "Xenium V1 Fresh Frozen Mouse Brain rep3 + ABCA Mouse Brain scRNAseq"
    dataset_description: "Xenium V1 Fresh Frozen Mouse Brain rep3 + ABCA Mouse Brain scRNAseq"
    dataset_organism: "mus_musculus"

  - id: "2023_10x_human_lung_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2023_10x_human_lung_xenium/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2020Travaglini_human_lung_sc/dataset.h5ad"
    dataset_id: "2023_10x_human_lung_xenium_combined"
    dataset_name: "Human lung combined 2023 10x Xenium 2020 Travaglini scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard"
    dataset_reference: "https://doi.org/10.1038/s41586-020-2922-4"
    dataset_summary: "Xenium Preview Human Non diseased Lung FFPE + 2020 Travaglini scRNAseq"
    dataset_description: "Xenium Preview Human Non diseased Lung FFPE + 2020 Travaglini scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2023_10x_human_lung_cancer_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2023_10x_human_lung_cancer_xenium/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "2023_10x_human_lung_cancer_xenium_combined"
    dataset_name: "Human lung cancer combined 2023 10x Xenium 2024 Zuani scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Xenium Preview Human Lung Cancer FFPE + 2024 Zuani scRNAseq"
    dataset_description: "Xenium Preview Human Lung Cancer FFPE + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2024_10x_human_skin_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2024_10x_human_skin_xenium/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2024Ganier_human_skin_sc/dataset.h5ad"
    dataset_id: "2024_10x_human_skin_xenium_combined"
    dataset_name: "Human skin combined 2024 10x Xenium 2024 Ganier scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/human-skin-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_reference: "https://doi.org/10.1073/pnas.2313326120"
    dataset_summary: "Xenium V1 hSkin nondiseased FFPE + 2024 Ganier scRNAseq"
    dataset_description: "Xenium V1 hSkin nondiseased FFPE + 2024 Ganier scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2024_10x_human_liver_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2024_10x_human_liver_xenium/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Andrews_human_liver_sc/dataset.h5ad"
    dataset_id: "2024_10x_human_liver_xenium_combined"
    dataset_name: "Human liver combined 2024 10x Xenium 2022 Andrews scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_reference: "https://doi.org/10.1002/hep4.1854"
    dataset_summary: "Xenium V1 hLiver FFPE + 2022 Andrews scRNAseq"
    dataset_description: "Xenium V1 hLiver FFPE + 2022 Andrews scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2024_10x_human_liver_cancer_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2024_10x_human_liver_cancer_xenium/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Lu_human_liver_cancer_sc/dataset.h5ad"
    dataset_id: "2024_10x_human_liver_cancer_xenium_combined"
    dataset_name: "Human liver cancer combined 2024 10x Xenium 2022 Lu scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_reference: "https://doi.org/10.1038/s41467-022-32283-3"
    dataset_summary: "Xenium V1 hLiver cancer FFPE + 2022 Lu scRNAseq"
    dataset_description: "Xenium V1 hLiver cancer FFPE + 2022 Lu scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2023_10x_human_colon_cancer_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2023_10x_human_colon_cancer_xenium/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2020Lee_human_colon_cancer_sc/dataset.h5ad"
    dataset_id: "2023_10x_human_colon_cancer_xenium_combined"
    dataset_name: "Human colon cancer combined 2023 10x Xenium 2020 Lee scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/human-colon-preview-data-xenium-human-colon-gene-expression-panel-1-standard"
    dataset_reference: "https://doi.org/10.1038/s41588-020-0636-z"
    dataset_summary: "Xenium V1 hColon Cancer FFPE + 2020 Lee scRNAseq"
    dataset_description: "Xenium V1 hColon Cancer FFPE + 2020 Lee scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2023_10x_human_breast_cancer_xenium_combined"
    input_sp: "$input_dir/10x_xenium/2023_10x_human_breast_cancer_xenium/dataset.zarr"
    input_sc: "$input_dir/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc/dataset.h5ad"
    dataset_id: "2023_10x_human_breast_cancer_xenium_combined"
    dataset_name: "Human breast cancer combined 2023 10x Xenium 2021 Wu scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard"
    dataset_reference: "https://doi.org/10.1038/s41588-021-00911-1"
    dataset_summary: "Xenium V1 FFPE Human Breast IDC + 2021 Wu scRNAseq"
    dataset_description: "Xenium V1 FFPE Human Breast IDC + 2021 Wu scRNAseq"
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
  --params-file /tmp/params_xenium.yaml \
  --config src/base/labels_nebius.config \
  --labels "task_ist_preprocessing,process_datasets,xenium"