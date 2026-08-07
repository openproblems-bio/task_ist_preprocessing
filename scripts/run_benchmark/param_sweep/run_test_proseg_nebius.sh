#!/bin/bash

# Nebius test run: all default methods + proseg at the transcript-assignment stage,
# with a parameter sweep over proseg's tuning knobs.
# See src/methods_transcript_assignment/proseg/NOTES.md ("Optimization / tuning").
#
# The stage default (basic_transcript_assignment) stays enabled alongside proseg so the run
# also produces the baseline the sweep is scored against — the workflow allows at most ONE
# non-default variant at a time, and proseg is that one non-default method here; every OTHER
# stage stays on its single default.
#
# PARAMS-FILE CAVEAT: `method_parameters_yaml` is opened by the WORKFLOW at runtime on the cloud
# (readYaml -> Nextflow file()), so a local /tmp path won't exist there and /scratch is read-only
# from the launch host. file() DOES stage http(s):// and this repo is public, so the sweep lives
# in a COMMITTED file read from GitHub via its raw URL => commit AND PUSH
# scripts/run_benchmark/param_sweep/proseg_params.yaml to $params_branch before launching.
# This is independent of --revision (which selects the pipeline CODE).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_test_s3="/scratch/task_ist_preprocessing/resources_test/task_ist_preprocessing/"
# Results publish to /scratch — created and written by the cloud compute env (read-only here).
publish_dir="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_proseg"

# The sweep lives in a committed file, read from GitHub at runtime. $params_branch defaults to
# the branch you are on; the file must be pushed there. (Independent of --revision below.)
params_repo="openproblems-bio/task_ist_preprocessing"
params_branch="$(git rev-parse --abbrev-ref HEAD)"
params_url="https://raw.githubusercontent.com/${params_repo}/${params_branch}/scripts/run_benchmark/param_sweep/proseg_params.yaml"

cat > /tmp/params_settings.yaml << HERE
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
transcript_assignment_methods:
  - basic_transcript_assignment
  - proseg
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
cat > /tmp/params.yaml << HERE
input_states: $resources_test_s3/**/state.yaml
rename_keys: 'input_sc:output_sc;input_sp:output_sp'
save_spatial_data: false
settings: '$(yq -o json /tmp/params_settings.yaml | jq -c .)'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# Fail early with a clear message if the params file isn't reachable on GitHub yet.
if ! curl -fsSL -o /dev/null "$params_url"; then
  echo "ERROR: params file not reachable at:" >&2
  echo "  $params_url" >&2
  echo "Commit and push scripts/run_benchmark/param_sweep/proseg_params.yaml to '$params_branch' first." >&2
  exit 1
fi

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/base/labels_nebius_test.config \
  --labels task_ist_preprocessing,test,proseg
