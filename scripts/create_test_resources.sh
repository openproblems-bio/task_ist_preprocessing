#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

###################################################################################
DATASET_ID="10x_xenium/2023_10x_mouse_brain_xenium"
TMP_DIR="temp/datasets/$DATASET_ID"
OUT_DIR="resources_test/common/2023_10x_mouse_brain_xenium"

# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip


rep1="$TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs"
rep2="$TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs"
rep3="$TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs"

if [ ! -d "$rep1" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip \
    -O $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
  unzip $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip -d $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_1
fi

if [ ! -d "$rep2" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_2/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip \
    -O $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip
  unzip $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs.zip -d $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_2
fi

if [ ! -d "$rep3" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_3/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip \
    -O $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip
  unzip $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs.zip -d $TMP_DIR/Xenium_V1_FF_Mouse_Brain_MultiSection_3
fi

# convert to zarr and concatenate
viash run src/data_loaders/download_10x_xenium/config.vsh.yaml -- \
  --input "$rep1" \
  --input "$rep2" \
  --input "$rep3" \
  --replicate_id rep1 \
  --replicate_id rep2 \
  --replicate_id rep3 \
  --output $TMP_DIR/full_dataset.zarr \
  --dataset_id "$DATASET_ID" \
  --dataset_name "Xenium V1 Fresh Frozen Mouse Brain" \
  --dataset_url "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard" \
  --dataset_summary "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1)." \
  --dataset_description "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1). Replicate results demonstrate the high reproducibility of data generated by the platform. 10x Genomics obtained tissue from a C57BL/6 mouse from Charles River Laboratories. Three adjacent 10µm sections were placed on the same slide. Tissues were prepared following the demonstrated protocols Xenium In Situ for Fresh Frozen Tissues - Tissue Preparation Guide (CG000579) and Xenium In Situ for Fresh Frozen Tissues - Fixation & Permeabilization (CG000581)." \
  --dataset_organism "mus_musculus" \
  --segmentation_id "cell;nucleus"

# crop the region
viash run src/data_processors/crop_region/config.vsh.yaml -- \
  --input "$TMP_DIR/full_dataset.zarr" \
  --output "$OUT_DIR/dataset.zarr" \
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


###################################################################################
DATASET_ID="allen_brain_cell_atlas/2023_Yao_mouse_brain_scRNAseq_10Xv2"
TMP_DIR="temp/datasets/$DATASET_ID"
OUT_DIR="resources_test/common/2023_abca_Yao_mouse_brain_scRNAseq_10Xv2"


# generate sc reference
VIASH_TEMP=/tmp/allen_brain_cell_atlas \
  viash run src/data_loaders/download_allen_brain_cell_atlas/config.vsh.yaml -- \
  --output "$TMP_DIR/tmp_sc_reference.h5ad" --regions "OLF;TH"

viash run src/data_processors/subset_reference/config.vsh.yaml -- \
  --input "$TMP_DIR/tmp_sc_reference.h5ad" \
  --output "$OUT_DIR/sc_reference.h5ad"



###################################################################################
aws s3 sync --profile op \
  "resources_test/common/2023_10x_mouse_brain_xenium" \
  "s3://openproblems-data/resources_test/common/2023_10x_mouse_brain_xenium" \
  --delete --dryrun
