#!/bin/bash

viash run src/core/fetch_task_execution/config.vsh.yaml -- \
  --input s3://openproblems-nextflow/work/f6/8565066aee4771cc2790b92b4ac660 \
  --output temp_debug_execution \
  --aws_profile op \
  --aws_credentials ~/.aws/credentials
