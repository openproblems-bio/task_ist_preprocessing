#!/bin/bash

# Build the combined Bruker CosMx lung-cancer *test* dataset used by
# results/Dataset_report.qmd.
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
# The SC input is the subsampled Zuani NSCLC test reference
# (scripts/create_test_resources/2024_zuani_human_nsclc_scrnaseq_nebius.sh).

set -e

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

SP=resources_test/common/2026_bruker_human_lung_cancer_cosmx/dataset.zarr
SC=resources_test/common/2024_zuani_human_nsclc_scrnaseq/dataset.h5ad
OUT_DIR=resources_test/task_ist_preprocessing/bruker_lung_cancer_cosmx_test_combined

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --dataset_id bruker_lung_cancer_cosmx_test_combined \
  --dataset_name "Test human lung cancer combined Bruker CosMx Lung9 rep1 2024 Zuani scRNAseq" \
  --dataset_url "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/" \
  --dataset_reference "" \
  --dataset_summary "Test crop of Bruker CosMx human lung cancer (Lung9 rep1) + subsampled 2024 Zuani NSCLC scRNAseq" \
  --dataset_description "Cropped test resource (Bruker CosMx human lung cancer, Lung9 rep1) combined with a subsampled Zuani 2024 NSCLC scRNAseq reference, for Dataset_report validation." \
  --dataset_organism "homo_sapiens" \
  --input_sc "$SC" \
  --input_sp "$SP" \
  --output_sc "$OUT_DIR/scrnaseq_reference.h5ad" \
  --output_sp "$OUT_DIR/raw_ist.zarr"

# The combined test dataset lands under $OUT_DIR (raw_ist.zarr has element "image").
# To publish it as a test resource, sync it to S3, e.g.:
#   aws s3 sync --profile op \
#     "$OUT_DIR" \
#     "s3://openproblems-data/resources_test/task_ist_preprocessing/bruker_lung_cancer_cosmx_test_combined" \
#     --delete --dryrun
