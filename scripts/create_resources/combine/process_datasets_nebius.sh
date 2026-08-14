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

# ── Batch 1: 10x Xenium ──────────────────────────────────────────────────────

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

  - id: "2026_10x_human_breast_cancer_atera_combined"
    input_sp: "$input_dir/10x_atera/2026_10x_human_breast_cancer_atera/dataset.zarr"
    input_sc: "$input_dir/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc/dataset.h5ad"
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

launch_batch /tmp/params_xenium.yaml "xenium"

# ── Batch 2: Vizgen MERSCOPE ──────────────────────────────────────────────────

cat > /tmp/params_vizgen.yaml << HERE
param_list:

  - id: "2022_vizgen_human_breast_cancer_merfish_combined/rep1"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_breast_cancer_merfish/rep1/dataset.zarr"
    input_sc: "$input_dir/wu_human_breast_cancer_sc/2021Wu_human_breast_cancer_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_breast_cancer_merfish_combined/rep1"
    dataset_name: "Human breast cancer combined 2022 Vizgen MERFISH rep1 2021 Wu scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://doi.org/10.1038/s41588-021-00911-1"
    dataset_summary: "Vizgen Human Breast Cancer MERFISH Patient1 + 2021 Wu scRNAseq"
    dataset_description: "Vizgen Human Breast Cancer MERFISH Patient1 + 2021 Wu scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_liver_cancer_merfish_combined/rep1"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_liver_cancer_merfish/rep1/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Lu_human_liver_cancer_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_liver_cancer_merfish_combined/rep1"
    dataset_name: "Human liver cancer combined 2022 Vizgen MERFISH rep1 2022 Lu scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://doi.org/10.1038/s41467-022-32283-3"
    dataset_summary: "Vizgen Human Liver Cancer MERFISH Patient1 + 2022 Lu scRNAseq"
    dataset_description: "Vizgen Human Liver Cancer MERFISH Patient1 + 2022 Lu scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_liver_cancer_merfish_combined/rep2"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_liver_cancer_merfish/rep2/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2022Lu_human_liver_cancer_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_liver_cancer_merfish_combined/rep2"
    dataset_name: "Human liver cancer combined 2022 Vizgen MERFISH rep2 2022 Lu scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://doi.org/10.1038/s41467-022-32283-3"
    dataset_summary: "Vizgen Human Liver Cancer MERFISH Patient2 + 2022 Lu scRNAseq"
    dataset_description: "Vizgen Human Liver Cancer MERFISH Patient2 + 2022 Lu scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_lung_cancer_merfish_combined/rep1"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_lung_cancer_merfish/rep1/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_lung_cancer_merfish_combined/rep1"
    dataset_name: "Human lung cancer combined 2022 Vizgen MERFISH rep1 2024 Zuani scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Vizgen Human Lung Cancer MERFISH Patient1 + 2024 Zuani scRNAseq"
    dataset_description: "Vizgen Human Lung Cancer MERFISH Patient1 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_lung_cancer_merfish_combined/rep2"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_lung_cancer_merfish/rep2/dataset.zarr"
    input_sc: "$input_dir/zuani_human_nsclc_sc/2024Zuani_human_nsclc_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_lung_cancer_merfish_combined/rep2"
    dataset_name: "Human lung cancer combined 2022 Vizgen MERFISH rep2 2024 Zuani scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526"
    dataset_summary: "Vizgen Human Lung Cancer MERFISH Patient2 + 2024 Zuani scRNAseq"
    dataset_description: "Vizgen Human Lung Cancer MERFISH Patient2 + 2024 Zuani scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_colon_cancer_merfish_combined/rep1"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_colon_cancer_merfish/rep1/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2020Lee_human_colon_cancer_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_colon_cancer_merfish_combined/rep1"
    dataset_name: "Human colon cancer combined 2022 Vizgen MERFISH rep1 2020 Lee scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://doi.org/10.1038/s41588-020-0636-z"
    dataset_summary: "2022 Vizgen Human Colon Cancer MERFISH Patient1 + 2020 Lee scRNAseq"
    dataset_description: "2022 Vizgen Human Colon Cancer MERFISH Patient1 + 2020 Lee scRNAseq"
    dataset_organism: "homo_sapiens"

  - id: "2022_vizgen_human_colon_cancer_merfish_combined/rep2"
    input_sp: "$input_dir/vizgen_merscope/2022_vizgen_human_colon_cancer_merfish/rep2/dataset.zarr"
    input_sc: "$input_dir/scrnaseq_for_ist/2020Lee_human_colon_cancer_sc/dataset.h5ad"
    dataset_id: "2022_vizgen_human_colon_cancer_merfish_combined/rep2"
    dataset_name: "Human colon cancer combined 2022 Vizgen MERFISH rep2 2020 Lee scRNAseq"
    dataset_url: "https://info.vizgen.com/ffpe-showcase"
    dataset_reference: "https://doi.org/10.1038/s41588-020-0636-z"
    dataset_summary: "2022 Vizgen Human Colon Cancer MERFISH Patient2 + 2020 Lee scRNAseq"
    dataset_description: "2022 Vizgen Human Colon Cancer MERFISH Patient2 + 2020 Lee scRNAseq"
    dataset_organism: "homo_sapiens"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

launch_batch /tmp/params_vizgen.yaml "vizgen"

# ── Batch 3: Bruker CosMx ─────────────────────────────────────────────────────

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

launch_batch /tmp/params_bruker.yaml "bruker"
