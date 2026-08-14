#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Nebius/Seqera (tw launch) variant of 2026_10x_human_breast_cancer_atera.sh.
# Publishes the cropped test resource to scratch instead of the local resources_test/.
# The atera zip URL is passed straight through (Nextflow stages it), so no wget step.
publish_dir="/scratch/task_ist_preprocessing/resources_test/common"

cat > /tmp/params_atera_test.yaml << HERE
param_list:
  - id: 2026_10x_human_breast_cancer_atera
    input: https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Breast_Cancer/WTA_Preview_FFPE_Breast_Cancer_outs.zip
    segmentation_id: [cell, nucleus]
    dataset_name: Atera WTA FFPE Human Breast Cancer
    dataset_url: https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer
    dataset_summary: Preview dataset showcasing the pre-commercial Atera Whole Transcriptome Assay (WTA) applied to FFPE human breast cancer tissue.
    dataset_description: "This human FFPE breast cancer data showcases results using the pre-commercial version of the Atera Whole Transcriptome Assay (WTA). A single 5µm FFPE section of breast cancer tissue (DCIS Grade 3, T1c N0 M0) was analyzed, yielding 170,057 detected cells with a median of 2,116 transcripts per cell."
    dataset_organism: homo_sapiens
    dataset_reference: NA
    # crop_region_* are in the "global" coordinate system, which for Xenium/Atera
    # is IMAGE PIXELS (0.2125 um/px), not microns: a 1000-unit box is ~212.5 um.
    # This window (px) recenters onto dense tumor tissue (~572 cells); the old
    # 5000-6000 box landed on a near-empty corner of the section (~13 cells).
    crop_region_min_x: 17000
    crop_region_max_x: 18000
    crop_region_min_y: 18000
    crop_region_max_y: 19000

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
  --params-file /tmp/params_atera_test.yaml \
  --config src/base/labels_nebius.config \
  --labels test_resources,atera

# The cropped dataset lands under $publish_dir/2026_10x_human_breast_cancer_atera/.
# To publish it as a test resource afterwards, sync it from scratch to S3, e.g.:
#   aws s3 sync --profile op \
#     "$publish_dir/2026_10x_human_breast_cancer_atera" \
#     "s3://openproblems-data/resources_test/common/2026_10x_human_breast_cancer_atera" \
#     --delete --dryrun
