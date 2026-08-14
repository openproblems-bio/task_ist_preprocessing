#!/bin/bash

# Convert the Kuppe kidney MERFISH (Vizgen MERSCOPE) samples into standardized SpatialData
# dataset.zarr via the generic process_vizgen_merscope workflow, publishing the loader output
# to the local /scratch raw folder (like the other *_nebius.sh spatial scripts).
#
# Inputs are the 2D (z3-only) region_0 copies mirrored to S3 by
# scripts/create_resources/spatial/upload_kuppe_kidney_merfish_2d.sh. The vizgen_merscope loader
# reads only z-plane 3 (spatialdata_io.merscope(z_layers=3)), so the z3-only copy is complete for
# this pipeline. Each region_0 already ships a Vizgen-native cell_boundaries.parquet
# (GeoParquet, ZIndex 0..6) which spatialdata_io reads directly (ZIndex==0 slice).
#
# NOTE: both samples have >65535 cell boundaries (cancer ~303k, healthy ~399k at ZIndex 0), so
# this exercises the loader's chunked-rasterization path — make sure the launch uses a
# vizgen_merscope image built AFTER that fix reached build/main (--revision build/main --pull-latest).
#
# The downstream combine step must read its spatial inputs from the same publish_dir
# (set input_dir="/scratch/task_ist_preprocessing/raw" there), pairing each spatial sample with the
# matching Kuppe SC reference (kuppe_kidney_{cancer,healthy}_sc, see process_kuppe_kidney_sc_nebius.sh).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# 2D region_0 copies produced by upload_kuppe_kidney_merfish_2d.sh
raw_dir="s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney/merfish"
# store the loader output locally, mirroring the process_datasets layout ($id/)
publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_kuppe_kidney_merfish.yaml << HERE
param_list:

  - id: "vizgen_merscope/kuppe_kidney_cancer_merfish/rep1"
    dataset_id: "vizgen_merscope/kuppe_kidney_cancer_merfish/rep1"
    input: "$raw_dir/cancer/region_0"
    dataset_name: "Kuppe Kidney Cancer MERFISH"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "Kuppe lab Clear Cell Kidney Carcinoma MERFISH (Vizgen MERSCOPE, 500-gene Pan-Cancer panel)."
    dataset_description: "Clear Cell Kidney Carcinoma FFPE sample profiled with MERFISH on the Vizgen MERSCOPE using the 500-gene Pan-Cancer panel (Kuppe lab). Cells segmented with Cellpose Cyto2 on the cell-boundary stain. Paired with the kuppe_kidney_cancer_sc snRNA-seq reference. Note: the Pan-Cancer panel does not separate NKT vs T cells or pericytes vs vSMC."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "vizgen_merscope/kuppe_kidney_healthy_merfish/rep1"
    dataset_id: "vizgen_merscope/kuppe_kidney_healthy_merfish/rep1"
    input: "$raw_dir/healthy/region_0"
    dataset_name: "Kuppe Kidney Healthy MERFISH"
    dataset_url: ""
    dataset_reference: ""
    dataset_summary: "Kuppe lab healthy kidney MERFISH (Vizgen MERSCOPE, 140-gene Spapros panel)."
    dataset_description: "Healthy kidney FFPE sample profiled with MERFISH on the Vizgen MERSCOPE using a 140-gene panel designed with Spapros from Kuppe et al 2021 snRNA-seq data. Cells segmented with Cellpose Cyto2 on the cell-boundary stain. Paired with the kuppe_kidney_healthy_sc snRNA-seq reference."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_vizgen_merscope/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_kuppe_kidney_merfish.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,kuppe_kidney_merfish
