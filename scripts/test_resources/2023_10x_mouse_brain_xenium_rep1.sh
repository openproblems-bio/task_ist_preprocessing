#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

if [ ! -d temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium ]; then
  mkdir -p temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium
fi
if [ ! -f temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip ]; then
  wget -O temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip \
    https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
fi

cat > /tmp/params.yaml << HERE
param_list:
  - id: 2023_10x_mouse_brain_xenium_rep1
    input: temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
    replicate_id: rep1
    segmentation_id:
      - cell
      - nucleus
    dataset_name: Xenium V1 Fresh Frozen Mouse Brain rep1
    dataset_url: https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard
    dataset_summary: Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform.
    dataset_description: Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1).
    dataset_organism: mus_musculus
    crop_region_min_x: 10000
    crop_region_max_x: 12000
    crop_region_min_y: 10000
    crop_region_max_y: 12000

publish_dir: resources_test/common
output: '\$id/dataset.h5ad'
output_state: '\$id/state.yaml'
HERE

# convert to zarr and concatenate
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_tenx_xenium/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml

aws s3 sync --profile op \
  "resources_test/common/2023_10x_mouse_brain_xenium" \
  "s3://openproblems-data/resources_test/common/2023_10x_mouse_brain_xenium" \
  --delete --dryrun
