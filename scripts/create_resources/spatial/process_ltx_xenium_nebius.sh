#!/bin/bash

# Process the custom LTX human-lung Xenium dataset (spatial side) — TWO samples.
# Reuses the generic 10x Xenium loader (process_tenx_xenium workflow) to convert
# each raw Xenium OUTPUT DIRECTORY (stored unzipped in the S3 raw_data bucket)
# into a standardized SpatialData dataset.zarr, published to the local /scratch
# raw folder (mirroring the layout process_datasets_*.sh reads from).
#
# Unlike the MPII run (which pointed at a single .zip), the LTX samples are stored
# as unzipped Xenium output directories, so `input` is the S3 directory prefix.
# The tenx_xenium loader handles a directory input: it walks the staged tree to
# find cell_feature_matrix.h5 and reads it with spatialdata_io.xenium().
#
# Each sample carries a multi-channel morphology_focus (ch0=DAPI); the loader keeps
# it as `morphology_mip`, which process_dataset later renames to the API-required
# `image`. Writing that multiscale image requires spatialdata>=0.8.0 in the loader
# container (older spatialdata + zarr>=3.2 fails with a rectilinear chunk-grid
# error). See the tenx_atera loader for the proven pin.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# store the loader output locally, mirroring the process_datasets layout ($id/)
# but under a sibling raw/ folder
publish_dir="/scratch/task_ist_preprocessing/raw"

# raw Xenium output directories in the S3 raw_data bucket (unzipped)
ltx_base="s3://openproblems-data/resources/raw_data/txSim_custom/ltx"

cat > /tmp/params_ltx_xenium.yaml << HERE
param_list:

  - id: "ltx_human_lung_xenium/B_AK-16-12103"
    input: "$ltx_base/B_AK-16-12103"
    dataset_name: "LTX human lung Xenium (sample B_AK-16-12103)"
    dataset_url: ""
    dataset_summary: "Custom LTX Xenium In Situ run of human lung tissue (sample B_AK-16-12103)."
    dataset_description: "Custom Xenium In Situ human lung dataset (LTX cohort, sample B_AK-16-12103). Paired with the LTX human lung scRNAseq annotation reference (ltx_full_annotation_2026_04_25)."
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "ltx_human_lung_xenium/L_AK14_14254"
    input: "$ltx_base/L_AK14_14254"
    dataset_name: "LTX human lung Xenium (sample L_AK14_14254)"
    dataset_url: ""
    dataset_summary: "Custom LTX Xenium In Situ run of human lung tissue (sample L_AK14_14254)."
    dataset_description: "Custom Xenium In Situ human lung dataset (LTX cohort, sample L_AK14_14254). Paired with the LTX human lung scRNAseq annotation reference (ltx_full_annotation_2026_04_25)."
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

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
  --params-file /tmp/params_ltx_xenium.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,ltx_xenium
