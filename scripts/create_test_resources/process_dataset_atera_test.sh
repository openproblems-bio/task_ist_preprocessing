#!/bin/bash

# Build the combined atera *test* dataset used by results/Dataset_report.qmd.
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
# The SC input is the subsampled Wu test reference (tw run 4DcTrWjb5pkW2u).

set -e

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

SP=resources_test/common/2026_10x_human_breast_cancer_atera/dataset.zarr
SC=resources_test/common/2021_wu_human_breast_cancer_scrnaseq/dataset.h5ad
OUT_DIR=resources_test/task_ist_preprocessing/atera_breast_cancer_test_combined

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --dataset_id atera_breast_cancer_test_combined \
  --dataset_name "Test atera WTA human breast cancer combined 2021 Wu scRNAseq" \
  --dataset_url "" \
  --dataset_reference "" \
  --dataset_summary "Test crop of Atera WTA FFPE human breast cancer + subsampled Wu 2021 breast scRNAseq" \
  --dataset_description "Cropped test resource (Atera WTA FFPE human breast cancer) combined with a subsampled Wu 2021 breast scRNAseq reference, for Dataset_report validation." \
  --dataset_organism "homo_sapiens" \
  --input_sc "$SC" \
  --input_sp "$SP" \
  --output_sc "$OUT_DIR/scrnaseq_reference.h5ad" \
  --output_sp "$OUT_DIR/raw_ist.zarr"

# The combined test dataset lands under $OUT_DIR (raw_ist.zarr has element "image").
# To publish it as a test resource, sync it to S3, e.g.:
#   aws s3 sync --profile op \
#     "$OUT_DIR" \
#     "s3://openproblems-data/resources_test/task_ist_preprocessing/atera_breast_cancer_test_combined" \
#     --delete --dryrun
