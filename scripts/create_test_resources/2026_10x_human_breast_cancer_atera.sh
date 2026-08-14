#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

if [ ! -d temp/datasets/10x_atera/2026_10x_human_breast_cancer_atera ]; then
  mkdir -p temp/datasets/10x_atera/2026_10x_human_breast_cancer_atera
fi
if [ ! -f temp/datasets/10x_atera/2026_10x_human_breast_cancer_atera/WTA_Preview_FFPE_Breast_Cancer_xe_outs.zip ]; then
  wget -O temp/datasets/10x_atera/2026_10x_human_breast_cancer_atera/WTA_Preview_FFPE_Breast_Cancer_xe_outs.zip \
     https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Breast_Cancer/WTA_Preview_FFPE_Breast_Cancer_outs.zip
fi

cat > /tmp/params.yaml << HERE
param_list:
  - id: 2026_10x_human_breast_cancer_atera
    input: temp/datasets/10x_atera/2026_10x_human_breast_cancer_atera/WTA_Preview_FFPE_Breast_Cancer_xe_outs.zip
    segmentation_id:
      - cell
      - nucleus
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

publish_dir: resources_test/common
output_dataset: '\$id/dataset.zarr'
output_state: '\$id/state.yaml'
HERE

# convert to zarr
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_tenx_atera/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml

# sync to s3
aws s3 sync --profile op \
  "resources_test/common/2026_10x_human_breast_cancer_atera" \
  "s3://openproblems-data/resources_test/common/2026_10x_human_breast_cancer_atera" \
  --delete --dryrun