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

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_vizgen.yaml \
  --config src/base/labels_nebius.config \
  --labels "task_ist_preprocessing,process_datasets,vizgen"