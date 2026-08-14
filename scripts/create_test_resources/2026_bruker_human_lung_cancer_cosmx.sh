#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Smallest of the CosMx archives mirrored to S3 (Lung9 Rep1: ~1.0 GiB flat + ~23.3 GiB
# morphology). The bruker_cosmx_nsclc loader streams these straight from S3 and only
# extracts what sopa needs, so no local download step is required (input_raw /
# input_morphology are plain strings, not staged files).
raw_dir="s3://openproblems-data/resources/raw_data/bruker_cosmx"

cat > /tmp/params.yaml << HERE
param_list:
  - id: 2026_bruker_human_lung_cancer_cosmx
    input_raw: $raw_dir/Lung9_Rep1+SMI+Flat+data.tar.gz
    input_morphology: $raw_dir/Lung9_Rep1+RawMorphologyImages.tar.gz
    segmentation_id:
      - cell
    dataset_name: Bruker CosMx Human Lung Cancer Lung9 Rep1
    dataset_url: https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
    dataset_summary: Bruker CosMx Human Lung Cancer Lung9 Rep1 dataset on FFPE.
    dataset_description: "Bruker CosMx Human Lung Cancer Lung9 Rep1 dataset on FFPE. Adenocarcinoma, G3, T3N1M0, IIIA, 65% tumour content."
    dataset_organism: homo_sapiens
    dataset_reference: NA
    # NOTE: coordinates are in the dataset's global coordinate system. Adjust the box
    # if it lands outside the tissue (an out-of-range crop yields an empty dataset).
    crop_region_min_x: 5000
    crop_region_max_x: 6000
    crop_region_min_y: 5000
    crop_region_max_y: 6000

publish_dir: resources_test/common
output_dataset: '\$id/dataset.zarr'
output_state: '\$id/state.yaml'
HERE

# convert to zarr
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_bruker_cosmx_nsclc/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml

# sync to s3
aws s3 sync --profile op \
  "resources_test/common/2026_bruker_human_lung_cancer_cosmx" \
  "s3://openproblems-data/resources_test/common/2026_bruker_human_lung_cancer_cosmx" \
  --delete --dryrun
