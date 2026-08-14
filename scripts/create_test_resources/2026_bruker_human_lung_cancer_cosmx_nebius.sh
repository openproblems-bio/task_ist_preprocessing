#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Nebius/Seqera (tw launch) variant of 2026_bruker_human_lung_cancer_cosmx.sh.
# Publishes the cropped test resource to scratch instead of the local resources_test/.
# Smallest CosMx archive pair (Lung9 Rep1); the bruker_cosmx_nsclc loader streams both
# straight from S3 and extracts only what sopa needs.
publish_dir="/scratch/task_ist_preprocessing/resources_test/common"
raw_dir="s3://openproblems-data/resources/raw_data/bruker_cosmx"

cat > /tmp/params_cosmx_test.yaml << HERE
param_list:
  - id: 2026_bruker_human_lung_cancer_cosmx
    input_raw: $raw_dir/Lung9_Rep1+SMI+Flat+data.tar.gz
    input_morphology: $raw_dir/Lung9_Rep1+RawMorphologyImages.tar.gz
    segmentation_id: ["cell"]
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

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_bruker_cosmx_nsclc/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_cosmx_test.yaml \
  --config src/base/labels_nebius.config \
  --labels test_resources,bruker_cosmx_nsclc

# The cropped dataset lands under $publish_dir/2026_bruker_human_lung_cancer_cosmx/.
# To publish it as a test resource afterwards, sync it from scratch to S3, e.g.:
#   aws s3 sync --profile op \
#     "$publish_dir/2026_bruker_human_lung_cancer_cosmx" \
#     "s3://openproblems-data/resources_test/common/2026_bruker_human_lung_cancer_cosmx" \
#     --delete --dryrun
