#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Nebius/Seqera (tw launch) variant of 2023_10x_mouse_brain_xenium_rep1.sh.
# Publishes the cropped test resource to scratch instead of the local resources_test/.
# The tenx_xenium loader accepts a directory, a zip, or a download url and stages/extracts
# it itself, so the 10x zip URL is passed straight through (no local wget/unzip step).
publish_dir="/scratch/task_ist_preprocessing/resources_test/common"

cat > /tmp/params_mouse_brain_xenium_test.yaml << HERE
param_list:
  - id: 2023_10x_mouse_brain_xenium_rep1
    input: https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
    segmentation_id: [cell, nucleus]
    dataset_name: Xenium V1 Fresh Frozen Mouse Brain rep1
    dataset_url: https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard
    dataset_summary: Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform.
    dataset_description: Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1).
    dataset_organism: mus_musculus
    dataset_reference: NA
    crop_region_min_x: 10000
    crop_region_max_x: 11000
    crop_region_min_y: 10000
    crop_region_max_y: 11000

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
  --params-file /tmp/params_mouse_brain_xenium_test.yaml \
  --config src/base/labels_nebius.config \
  --labels test_resources,tenx_xenium

# The cropped dataset lands under $publish_dir/2023_10x_mouse_brain_xenium_rep1/.
# To publish it as a test resource afterwards, sync it from scratch to S3, e.g.:
#   aws s3 sync --profile op \
#     "$publish_dir/2023_10x_mouse_brain_xenium_rep1" \
#     "s3://openproblems-data/resources_test/common/2023_10x_mouse_brain_xenium_rep1" \
#     --delete --dryrun
