#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_OUT=resources/datasets_raw/10x_fresh_frozen_mouse_brain_replicates
RESOURCES_OUT=resources/datasets/10x_xenium/10x_fresh_frozen_mouse_brain_replicates

# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip


rep1="$RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs"
if [ ! -d "$rep1" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip \
    -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
  unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1
fi

rep2="$RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs"
if [ ! -d "$rep2" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip \
    -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
  unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2
fi

rep3="$RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs"
if [ ! -d "$rep3" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip \
    -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip
  unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3
fi

# convert to zarr and concatenate
viash run src/data_loaders/download_10x_xenium/config.vsh.yaml -- \
  --input "$rep1" \
  --input "$rep2" \
  --input "$rep3" \
  --replicate_id rep1 \
  --replicate_id rep2 \
  --replicate_id rep3 \
  --output $RAW_OUT/full_dataset.zarr

# crop the region
viash run src/data_processors/crop_region/config.vsh.yaml -- \
  --input $RAW_OUT/full_dataset.zarr \
  --output $RESOURCES_OUT/dataset.zarr \
  --replicate_id "rep1" \
  --min_x 10000 \
  --max_x 12000 \
  --min_y 10000 \
  --max_y 12000 \
  --replicate_id rep2 \
  --min_x 10000 \
  --max_x 12000 \
  --min_y 10000 \
  --max_y 12000 \
  --replicate_id rep3 \
  --min_x 10000 \
  --max_x 12000 \
  --min_y 10000 \
  --max_y 12000
