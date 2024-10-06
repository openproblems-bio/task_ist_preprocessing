#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

cat > /tmp/params.yaml << HERE
param_list:
  - id: 2023_yao_mouse_brain_scrnaseq_10xv2
    regions:
      - OLF
      - TH
    dataset_name: ABCA Mouse Brain scRNAseq
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717
    dataset_reference: 10.1038/s41586-023-06812-z
    dataset_summary: A high-resolution scRNAseq atlas of cell types in the whole mouse brain
    dataset_description: See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset.
    dataset_organism: mus_musculus

do_subsample: true
n_obs: 400
n_vars: 10000

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: resources_test/common
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_allen_brain_cell_atlas/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml

aws s3 sync --profile op \
  "resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2" \
  "s3://openproblems-data/resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2" \
  --delete --dryrun
