#!/bin/bash

# Process ONLY the MPII custom human-lung Xenium dataset (spatial side).
# Converts the raw Xenium output zip (stored in the S3 raw_data bucket) into the
# standardized SpatialData dataset.zarr via the process_tenx_xenium workflow,
# and publishes it to the local /scratch raw folder (mirroring the layout that
# process_datasets_mpii_nebius.sh reads from).
#
# NOTE: this is a genuine Xenium run (instrument XETG00117). If the run fails
# with an irregular / rectilinear chunk-grid error while reading the loader
# output, switch the --main-script below to the tenx_atera workflow
# (target/nextflow/datasets/workflows/process_tenx_atera/main.nf), whose loader
# carries the rechunk_uniform() fix for exactly that issue.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# store the loader output locally, mirroring the process_datasets layout ($id/)
# but under a sibling raw/ folder
publish_dir="/scratch/task_ist_preprocessing/raw"

# raw Xenium output zip in the S3 raw_data bucket
input_zip="s3://openproblems-data/resources/raw_data/txSim_custom/MPII/978_reg1_output-XETG00117__0015978__Region_1__20240718__175145.zip"

cat > /tmp/params_mpii_xenium.yaml << HERE
param_list:

  - id: "mpii_human_lung_xenium/978_reg1"
    input: "$input_zip"
    dataset_name: "MPII human lung Xenium (sample 978, region 1)"
    dataset_url: ""
    dataset_summary: "Custom MPII Xenium In Situ run of human lung tissue (sample 978, region 1, XETG00117, 2024-07-18)."
    dataset_description: "Custom Xenium In Situ human lung dataset generated at MPII (sample 978, region 1, instrument XETG00117, acquired 2024-07-18). Paired with the MPII human lung scRNAseq annotation reference."
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_tenx_xenium/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_mpii_xenium.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,mpii_xenium
