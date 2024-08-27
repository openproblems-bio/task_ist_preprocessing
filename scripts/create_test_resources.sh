#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_OUT=resources/datasets_raw/10x_fresh_frozen_mouse_brain_replicates

# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip


wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip \
  -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_1

wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip \
  -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_2

wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip \
  -O $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip
unzip $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip -d $RAW_OUT/Xenium_V1_FF_Mouse_Brain_MultiSection_3
