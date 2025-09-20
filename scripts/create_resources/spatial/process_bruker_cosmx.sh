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
    input_raw: "https://smi-public.objects.liquidweb.services/HalfBrain.zip"
    input_flat_files: "https://smi-public.objects.liquidweb.services/Half%20%20Brain%20simple%20%20files%20.zip"
    dataset_name: "Bruker CosMx Mouse Brain"
    dataset_url: "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-smi-mouse-brain-ffpe-dataset/"
    dataset_summary: "Bruker CosMx Mouse Brain dataset on FFPE covering a full hemisphere of a mouse brain."
    dataset_description: "Bruker CosMx Mouse Brain dataset on FFPE covering a full hemisphere of a mouse brain."
    dataset_organism: "mus_musculus"
    segmentation_id: ["cell"]


output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_bruker_cosmx/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,bruker_cosmx