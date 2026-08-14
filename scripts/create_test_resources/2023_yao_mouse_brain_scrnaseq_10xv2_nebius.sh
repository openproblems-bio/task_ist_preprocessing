#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# Nebius/Seqera (tw launch) variant of 2023_yao_mouse_brain_scrnaseq_10xv2.sh.
# Matching single-cell reference for the 10x Xenium mouse-brain spatial test sample.
# Publishes to scratch instead of the local resources_test/.
#
# The allen_brain_cell_atlas loader downloads the requested regions from the ABCA directly
# (no local input), and subsamples during loading. This workflow does NOT use the
# processors/subsample component, so it is a single step (unlike the Wu/Zuani SC chains).
#
# NOTE: the current workflow subsamples via the loader's --sample_n_obs (cells only); the
# older do_subsample/n_obs/n_vars params in the local script are no longer part of the
# workflow interface and are silently ignored by findStates (which would yield the FULL
# OLF+TH regions rather than a small test reference). Use sample_n_obs here.
publish_dir="/scratch/task_ist_preprocessing/resources_test/common"

cat > /tmp/params_yao_mouse_brain_test.yaml << HERE
param_list:
  - id: 2023_yao_mouse_brain_scrnaseq_10xv2
    regions: [OLF, TH]
    dataset_name: ABCA Mouse Brain scRNAseq
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717
    dataset_reference: 10.1038/s41586-023-06812-z
    dataset_summary: A high-resolution scRNAseq atlas of cell types in the whole mouse brain
    dataset_description: See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset.
    dataset_organism: mus_musculus
    sample_n_obs: 400
    sample_seed: 0

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
  --params-file /tmp/params_yao_mouse_brain_test.yaml \
  --config src/base/labels_nebius.config \
  --labels test_resources,allen_brain_cell_atlas

# The subsampled reference lands under $publish_dir/2023_yao_mouse_brain_scrnaseq_10xv2/.
# To publish it as a test resource afterwards, sync it from scratch to S3, e.g.:
#   aws s3 sync --profile op \
#     "$publish_dir/2023_yao_mouse_brain_scrnaseq_10xv2" \
#     "s3://openproblems-data/resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2" \
#     --delete --dryrun
