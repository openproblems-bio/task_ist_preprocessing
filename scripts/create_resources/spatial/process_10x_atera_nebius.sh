#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# store the loader output locally, mirroring the process_datasets layout ($id/)
# but under a sibling raw/ folder
publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_atera.yaml << HERE
param_list:

  - id: "10x_atera/2026_10x_human_breast_cancer_atera"
    input:  https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Breast_Cancer/WTA_Preview_FFPE_Breast_Cancer_outs.zip
    dataset_name: "Atera WTA FFPE Human Breast Cancer"
    dataset_url: "https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer"
    dataset_summary: "Preview dataset showcasing the pre-commercial Atera Whole Transcriptome Assay (WTA) applied to FFPE human breast cancer tissue, profiling 18,028 genes and detecting 170,057 cells."
    dataset_description: "This human FFPE breast cancer data showcases results using the pre-commercial version of the Atera Whole Transcriptome Assay (WTA), which is currently under development. The assay is designed to closely match the Chromium Flex Apex assay in terms of content and sensitivity, and includes 18,028 genes. A single 5µm FFPE section of breast cancer tissue (DCIS Grade 3, T1c N0 M0) was analyzed, yielding 170,057 detected cells with a median of 2,116 transcripts per cell and 624,095,990 total high-quality decoded transcripts across 58.9 million µm² of tissue area. Output files are formatted to closely resemble Xenium Onboard Analysis v4 file formats."
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_tenx_atera/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_atera.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,atera