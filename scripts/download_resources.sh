#!/bin/bash

set -e

echo ">> Downloading resources"

# the sync_resources script uses the test_resources S3 URI's in the _viash.yaml to download the resources.
common/scripts/sync_resources --delete
