#!/bin/bash

# Nebius test run: all default methods + StarDist2D (stardist) segmentation, with a
# parameter sweep over the optimization levers identified for StarDist.
# See src/methods_segmentation/stardist/NOTES.md ("Optimization / tuning").
# Local sibling: run_test_stardist_local.sh
#
# stardist runs on TensorFlow and is scheduled with the `gpu` label (its config uses
# a GPU-capable base image), hence the Nebius GPU compute env.
#
# PARAMS-FILE CAVEAT (why this differs from the local script):
#   `tw launch --params-file` is read client-side, but `method_parameters_yaml`
#   is a path the WORKFLOW opens at runtime on the cloud (readYaml -> Nextflow
#   file()). A local /tmp path does not exist there, and /scratch (where results
#   publish) is READ-ONLY from the launch host — which is why the binning
#   method_params block is commented out in run_test_nebius.sh. file() does stage
#   http(s):// though, and this repo is public, so we keep the sweep in a COMMITTED
#   file (scripts/run_benchmark/stardist_params.yaml) and read it from GitHub via
#   its raw URL. => the params file must be committed AND PUSHED to $params_branch
#   before launching (edit the file there, not here, to change the sweep).

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_test_s3=s3://openproblems-data/resources_test/task_ist_preprocessing
# Results publish to /scratch — created and written by the cloud compute env, so
# the launcher does NOT create it here (it is read-only from the launch host).
publish_dir="/scratch/results/runs/$(date +%Y-%m-%d_%H-%M-%S)_stardist"

# The sweep lives in a committed file, read from GitHub at runtime. $params_branch
# defaults to the branch you are on; the file must be pushed there on GitHub. (This
# is independent of --revision below, which selects the pipeline CODE to run.)
params_repo="openproblems-bio/task_ist_preprocessing"
params_branch="$(git rev-parse --abbrev-ref HEAD)"
params_url="https://raw.githubusercontent.com/${params_repo}/${params_branch}/scripts/run_benchmark/stardist_params.yaml"

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
  - stardist
#  - cellpose
#  - cellposev4
#  - binning
#  - watershed
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
  echo "Commit and push scripts/run_benchmark/stardist_params.yaml to '$params_branch' first." >&2
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
  --config src/base/labels_nebius.config \
  --labels task_ist_preprocessing,test,stardist,gpu
