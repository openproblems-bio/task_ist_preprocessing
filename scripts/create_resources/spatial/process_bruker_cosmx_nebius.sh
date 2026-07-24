#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "bruker_cosmx/bruker_mouse_brain_cosmx/rep1"
    input_raw: "s3://openproblems-data/resources/raw_data/bruker_cosmx/HalfBrain.zip"
    input_flat_files: "s3://openproblems-data/resources/raw_data/bruker_cosmx/Half Brain simple files.zip"
    dataset_name: "Bruker CosMx Mouse Brain"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-smi-mouse-brain-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Mouse Brain dataset on FFPE covering a full hemisphere of a mouse brain."
    dataset_description: "Bruker CosMx Mouse Brain dataset on FFPE covering a full hemisphere of a mouse brain."
    dataset_organism: "mus_musculus"
    segmentation_id: ["cell"]

  - id: "bruker_cosmx/bruker_human_liver_cosmx"
    input_raw: "s3://openproblems-data/resources/raw_data/bruker_cosmx/NormalLiverFiles.zip"
    dataset_name: "Bruker CosMx Human Liver"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Human Liver dataset on FFPE."
    dataset_description: "Bruker CosMx Human Liver dataset on FFPE."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]


output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_bruker_cosmx/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,bruker_cosmx
