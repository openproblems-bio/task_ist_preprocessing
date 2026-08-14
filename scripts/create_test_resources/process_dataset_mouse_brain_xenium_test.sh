#!/bin/bash

# Build the combined mouse-brain Xenium *test* dataset used by results/Dataset_report.qmd.
#
# Why this exists: the loaders emit the image element as `morphology_mip`; only
# process_dataset renames it to the API-required `image` (file_common_ist.yaml).
# The resources_test loader outputs (resources_test/common/<id>/dataset.zarr) are
# raw loader outputs, so they keep `morphology_mip` and lack a panel-subset SC —
# the report can't read them directly. Running process_dataset here is the
# workaround: it renames morphology_mip -> image, subsets the SC to the shared
# panel, and writes a combined dataset (raw_ist.zarr + scrnaseq_reference.h5ad),
# exactly like test_pipeline.sh does for mouse_brain_combined.
#
# The SC input is the subsampled Yao mouse-brain test reference
# (scripts/create_test_resources/2023_yao_mouse_brain_scrnaseq_10xv2_nebius.sh).

set -e

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

SP=resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr
SC=resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad
OUT_DIR=resources_test/task_ist_preprocessing/mouse_brain_xenium_test_combined

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --dataset_id mouse_brain_xenium_test_combined \
  --dataset_name "Test mouse brain combined 2023 10x Xenium rep1 2023 Yao scRNAseq" \
  --dataset_url "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard" \
  --dataset_reference "" \
  --dataset_summary "Test crop of 10x Xenium fresh-frozen mouse brain (rep1) + subsampled 2023 Yao mouse brain scRNAseq" \
  --dataset_description "Cropped test resource (10x Xenium fresh-frozen mouse brain, rep1) combined with a subsampled Yao 2023 mouse brain scRNAseq reference, for Dataset_report validation." \
  --dataset_organism "mus_musculus" \
  --input_sc "$SC" \
  --input_sp "$SP" \
  --output_sc "$OUT_DIR/scrnaseq_reference.h5ad" \
  --output_sp "$OUT_DIR/raw_ist.zarr"

# The combined test dataset lands under $OUT_DIR (raw_ist.zarr has element "image").
# To publish it as a test resource, sync it to S3, e.g.:
#   aws s3 sync --profile op \
#     "$OUT_DIR" \
#     "s3://openproblems-data/resources_test/task_ist_preprocessing/mouse_brain_xenium_test_combined" \
#     --delete --dryrun
