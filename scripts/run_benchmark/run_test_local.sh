#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

echo "Running benchmark on test data"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="testrun_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="temp/results/${RUN_ID}"

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c common/nextflow_helpers/labels_ci.config \
  --id mouse_brain_combined \
  --input_sc resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad \
  --input_sp resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr \
  --output_state state.yaml \
  --publish_dir "$publish_dir"
