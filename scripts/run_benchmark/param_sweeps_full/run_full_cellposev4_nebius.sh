#!/bin/bash

# Nebius FULL run: all default methods + Cellpose 4 (cellposev4) segmentation, driven by
# the FOLLOW-UP (narrowed, Pareto-selected) parameter sweep, on the FULL datasets and the
# FULL (non-test) labels_nebius.config.
#
# Sweep file: param_sweeps_full/cellposev4_params.yaml  (1 default + 10 swept = 11 variants;
#   AUTO-GENERATED from the test-resource sweep results — do not hand-edit it here).
# Small/test sibling: ../param_sweep/run_test_cellposev4_nebius.sh
#
# cellposev4 (cpsam / SAM ViT) is GPU-heavy, hence the `gpu` label; GPU routing is by
# src/base/labels_nebius.config nodeSelectors on the component's `gpu` label.
#
# PARAMS-FILE CAVEAT: `method_parameters_yaml` is opened by the WORKFLOW at runtime on the
# cloud (readYaml -> Nextflow file()), which stages http(s):// but NOT a launch-host local
# path. So the sweep is read from GitHub via its raw URL => it must be committed AND PUSHED
# to $params_branch before launching. (Independent of --revision, which selects pipeline CODE.)

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# FULL datasets on the Nebius shared /scratch mount (not the small S3 test set).
resources_s3=/scratch/task_ist_preprocessing/datasets
publish_dir="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_cellposev4_full_sweep"

# The sweep lives in a committed file, read from GitHub at runtime. $params_branch
# defaults to the branch you are on; the file must be pushed there on GitHub.
params_repo="openproblems-bio/task_ist_preprocessing"
params_branch="$(git rev-parse --abbrev-ref HEAD)"
params_url="https://raw.githubusercontent.com/${params_repo}/${params_branch}/scripts/run_benchmark/param_sweeps_full/cellposev4_params.yaml"

cat > /tmp/params_settings_cellposev4_full.yaml << HERE
default_methods:
  - custom_segmentation
  - basic_transcript_assignment
  - basic_count_aggregation
  - basic_qc_filter
  - alpha_shapes
  - normalize_by_volume
  - tacco
  - no_correction
segmentation_methods:
  - custom_segmentation
  - cellposev4
transcript_assignment_methods:
  - basic_transcript_assignment
count_aggregation_methods:
  - basic_count_aggregation
qc_filtering_methods:
  - basic_qc_filter
volume_calculation_methods:
  - alpha_shapes
normalization_methods:
  - normalize_by_volume
celltype_annotation_methods:
  - tacco
expression_correction_methods:
  - no_correction
gene_efficiency_correction_methods:
  - no_correction
method_parameters_yaml: $params_url
HERE

# Write the parameters to file (input_states version, NOTE: enable `-entry_name auto` for this)
cat > /tmp/params_cellposev4_full.yaml << HERE
input_states: $resources_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings_cellposev4_full.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# Fail early with a clear message if the params file isn't reachable on GitHub yet.
if ! curl -fsSL -o /dev/null "$params_url"; then
  echo "ERROR: params file not reachable at:" >&2
  echo "  $params_url" >&2
  echo "Commit and push scripts/run_benchmark/param_sweeps_full/cellposev4_params.yaml to '$params_branch' first." >&2
  exit 1
fi

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params_cellposev4_full.yaml \
  --entry-name auto \
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,full,cellposev4,gpu
