#!/bin/bash

# Nebius version of process_allen_brain_cell_atlas_brain.sh: process the ABCA (2023 Yao)
# mouse-brain scRNAseq reference on Nebius and publish to the scratch raw/ folder (like the
# other *_nebius.sh scripts), so the scratch-reading combine scripts (e.g.
# process_datasets_bruker_nebius.sh) find it at
#   /scratch/task_ist_preprocessing/raw/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad
# The S3 sibling (process_allen_brain_cell_atlas_brain.sh) publishes to resources/datasets.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="/scratch/task_ist_preprocessing/raw"

cat > /tmp/params_allen_brain_sc.yaml << HERE
param_list:

  - id: allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2
    regions:
      - CTXsp
      - HPF
      - HY
      - Isocortex
      - MB
      - OLF
      - TH
    dataset_name: ABCA Mouse Brain scRNAseq
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717
    dataset_reference: 10.1038/s41586-023-06812-z
    dataset_summary: A high-resolution scRNAseq atlas of cell types in the whole mouse brain
    dataset_description: See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset.
    dataset_organism: mus_musculus

sample_n_obs: 500000
sample_obs_weight: subclass
sample_transform: log
sample_seed: 42
keep_files: false # disk isn't large enough

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_allen_brain_cell_atlas/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_allen_brain_sc.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,allen_brain_cell_atlas
