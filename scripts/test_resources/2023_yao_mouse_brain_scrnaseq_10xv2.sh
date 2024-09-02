#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASET_ID="2023_yao_mouse_brain_scrnaseq_10xv2"
TMP_DIR="temp/datasets/$DATASET_ID"
OUT_DIR="resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2"

# generate sc reference
VIASH_TEMP=/tmp/allen_brain_cell_atlas \
  viash run src/data_loaders/download_allen_brain_cell_atlas/config.vsh.yaml \
  --keep true -- \
  --regions "OLF;TH" \
  --output "$TMP_DIR/tmp_dataset.h5ad" \
  --dataset_id "$DATASET_ID" \
  --dataset_name "ABCA Mouse Brain scRNAseq" \
  --dataset_url "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717" \
  --dataset_reference "10.1038/s41586-023-06812-z" \
  --dataset_summary "A high-resolution scRNAseq atlas of cell types in the whole mouse brain" \
  --dataset_description "See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset." \
  --dataset_organism "mus_musculus"

viash run src/data_processors/subset_reference/config.vsh.yaml -- \
  --input "$TMP_DIR/tmp_dataset.h5ad" \
  --output "$OUT_DIR/dataset.h5ad"

aws s3 sync --profile op \
  "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2" \
  "s3://openproblems-data/resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2" \
  --delete --dryrun
