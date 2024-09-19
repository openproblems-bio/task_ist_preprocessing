#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resource/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2
    regions:
      - CTXsp
      - HPF
      - HY
      - Isocortex-1
      - Isocortex-2
      - Isocortex-3
      - Isocortex-4
      - MB
      - OLF
      - TH
    dataset_name: ABCA Mouse Brain scRNAseq
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717
    dataset_reference: 10.1038/s41586-023-06812-z
    dataset_summary: A high-resolution scRNAseq atlas of cell types in the whole mouse brain
    dataset_description: See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset.
    dataset_organism: mus_musculus

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch openproblems-bio/task_ist_preprocessing \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_allen_brain_cell_atlas/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,allen_brain_cell_atlas