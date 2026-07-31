#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

# Raw NSCLC archives mirrored to S3 (see scripts/create_resources/spatial/mirror_bruker_to_s3.sh).
# Each sample ships two tar.gz archives: the flat files + cell labels, and the raw morphology images.
raw_dir="s3://openproblems-data/resources/raw_data/bruker_cosmx"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep1"
    input_raw: "$raw_dir/Lung5_Rep1+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung5_Rep1+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung5 Rep1"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep1 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep1 dataset on FFPE. Adenocarcinoma, G1, T2aN2M0, IIIA, 75% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep2"
    input_raw: "$raw_dir/Lung5_Rep2+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung5_Rep2+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung5 Rep2"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep2 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep2 dataset on FFPE. Adenocarcinoma, G1, T2aN2M0, IIIA, 75% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung5_rep3"
    input_raw: "$raw_dir/Lung5_Rep3+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung5_Rep3+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung5 Rep3"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung5 Rep3 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung5 Rep3 dataset on FFPE. Adenocarcinoma, G1, T2aN2M0, IIIA, 75% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung6"
    input_raw: "$raw_dir/Lung6+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung6+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung6"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung6 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung6 dataset on FFPE. Squamous cell carcinoma, G2, T2bN2M0, IIIA, 90% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung9_rep1"
    input_raw: "$raw_dir/Lung9_Rep1+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung9_Rep1+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung9 Rep1"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung9 Rep1 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung9 Rep1 dataset on FFPE. Adenocarcinoma, G3, T3N1M0, IIIA, 65% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung9_rep2"
    input_raw: "$raw_dir/Lung9_Rep2+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung9_Rep2+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung9 Rep2"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung9 Rep2 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung9 Rep2 dataset on FFPE. Adenocarcinoma, G3, T3N1M0, IIIA, 65% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung12"
    input_raw: "$raw_dir/Lung12+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung12+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung12"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung12 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung12 dataset on FFPE. Adenocarcinoma, G3, T4N0M0, IIIA, 85% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_lung_cancer_cosmx/lung13"
    input_raw: "$raw_dir/Lung13+SMI+Flat+data.tar.gz"
    input_morphology: "$raw_dir/Lung13+RawMorphologyImages.tar.gz"
    dataset_name: "Bruker CosMx Human Lung Cancer Lung13"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Lung Cancer Lung13 dataset on FFPE."
    dataset_description: "Bruker CosMx Human Lung Cancer Lung13 dataset on FFPE. Adenocarcinoma, G1, T3N0M0, IIB, 55% tumour content."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]


output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_bruker_cosmx_nsclc/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,bruker_cosmx_nsclc
