#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~ USAGE ~~~~~~~~~~~~~~~~~~~~~~~#
# bash ./scripts/dataset_downloads/2023Janesick_human_breast_cancer_xenium.sh [OUTPUT_DIR]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# If an argument is passed, use it as the output directory, otherwise use the default
if [ -z "$1" ]; then
  OUT_DIR="resources_test/common/2023Janesick_human_breast_cancer_xenium"
else
  OUT_DIR="$1"
fi


# Download cell type annotations if not already present
if [ ! -f "$OUT_DIR/Cell_Barcode_Type_Matrices.xlsx" ]; then
  wget https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx \
    -O "$OUT_DIR/Cell_Barcode_Type_Matrices.xlsx"
fi

# Download spatial data zips and unzip
# https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip
# https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip

rep1="$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs"
rep2="$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs"

if [ ! -d "$rep1" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip \
    -O "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip"
  unzip "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip" -d "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep1"
fi

if [ ! -d "$rep2" ]; then
  wget https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip \
    -O "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip"
  unzip "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip" -d "$OUT_DIR/Xenium_FFPE_Human_Breast_Cancer_Rep2"
fi
    