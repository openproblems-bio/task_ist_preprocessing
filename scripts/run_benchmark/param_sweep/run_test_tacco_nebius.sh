#!/bin/bash

# Nebius test run: all default methods + TACCO (tacco) cell-type annotation, with a
# parameter sweep over the tacco.tl.annotate optimization knobs.
# See src/methods_cell_type_annotation/tacco/NOTES.md ("Optimization / tuning").
#
# TACCO (optimal-transport annotation transfer) is CPU-only — NO GPU, so there is no
# `gpu` label here (this mirrors run_test_nebius.sh / the singler launcher). The compute
# env is the same one all these test runs use.
#
# The annotation stage keeps its default method `tacco` alongside the swept `tacco`
# variants so the run has a baseline; every other stage is left on its single default
# (the swept method must be the only non-default thing). Note tacco IS the stage default,
# so here the "baseline" and the swept method are the same component — the default
# variant reproduces stock tacco.
#
# !!! REQUIRES A REBUILD BEFORE THIS DOES ANYTHING !!!
#   The sweep varies NEWLY EXPOSED args (method, multi_center, platform_iterations,
#   bisections) that do NOT exist in the deployed build/main tacco container yet. Merge
#   the tacco config.vsh.yaml + script.py change to main, let CI rebuild the container,
#   and verify with `check-component tacco` BEFORE launching — otherwise every non-default
#   variant errors with "unknown option". run_sweep's build gate enforces this.
#
# PARAMS-FILE CAVEAT (why this reads from GitHub):
#   `tw launch --params-file` is read client-side, but `method_parameters_yaml` is a path
#   the WORKFLOW opens at runtime on the cloud (readYaml -> Nextflow file()). A local /tmp
#   path does not exist there, and /scratch (where results publish) is READ-ONLY from the
#   launch host. file() does stage http(s):// though, and this repo is public, so we keep
#   the sweep in a COMMITTED file (scripts/run_benchmark/param_sweep/tacco_params.yaml)
#   and read it from GitHub via its raw URL. => the params file must be committed AND
#   PUSHED to $params_branch before launching. This is independent of --revision below,
#   which selects the pipeline CODE.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_test_s3="/scratch/task_ist_preprocessing/resources_test/task_ist_preprocessing/"
# Results publish to /scratch — created and written by the cloud compute env, so
# the launcher does NOT create it here (it is read-only from the launch host).
publish_dir="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_tacco"

# The sweep lives in a committed file, read from GitHub at runtime. $params_branch
# defaults to the branch you are on; the file must be pushed there on GitHub. (This
# is independent of --revision below, which selects the pipeline CODE to run.)
params_repo="openproblems-bio/task_ist_preprocessing"
params_branch="$(git rev-parse --abbrev-ref HEAD)"
params_url="https://raw.githubusercontent.com/${params_repo}/${params_branch}/scripts/run_benchmark/param_sweep/tacco_params.yaml"

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
  echo "Commit and push scripts/run_benchmark/param_sweep/tacco_params.yaml to '$params_branch' first." >&2
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
  --labels task_ist_preprocessing,test,tacco
