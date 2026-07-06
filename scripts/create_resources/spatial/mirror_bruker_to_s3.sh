#!/bin/bash

# Mirror the Bruker CosMx raw zips from the (flaky) NanoString/liquidweb server to S3.
#
# The NanoString server frequently resets the connection mid-transfer, and these
# zips are enormous (HalfBrain ~175 GB, NormalLiver ~349 GB). Nextflow's foreign-file
# staging does NOT resume, so it dies on every run. Mirroring to S3 once makes future
# runs fast and deterministic.
#
# This script is resumable and idempotent:
#   - curl -C - resumes a partial download after a connection reset
#   - files already present in S3 (with matching size) are skipped
#   - each file is deleted from local scratch after upload, so you only ever need
#     free disk for the LARGEST single file (~350 GB), not the sum
#
# Just re-run it if it dies; it picks up where it left off.
#
# Usage:
#   SCRATCH_DIR=/path/to/big/scratch ./mirror_bruker_to_s3.sh

set -euo pipefail

# --- Config -----------------------------------------------------------------

SRC_BASE="https://smi-public.objects.liquidweb.services"
S3_DEST="${S3_DEST:-s3://openproblems-data/resources_raw/bruker_cosmx}"
SCRATCH_DIR="${SCRATCH_DIR:-$PWD/bruker_mirror_scratch}"

# Files to mirror. The URL-encoded names are what the server serves; the second
# column is the (decoded) name to store under on S3.
FILES=(
  "HalfBrain.zip|HalfBrain.zip"
  "Half%20%20Brain%20simple%20%20files%20.zip|Half Brain simple files.zip"
  "NormalLiverFiles.zip|NormalLiverFiles.zip"
)

# --- Preflight --------------------------------------------------------------

command -v curl >/dev/null || { echo "ERROR: curl not found" >&2; exit 1; }
command -v aws  >/dev/null || { echo "ERROR: aws CLI not found" >&2; exit 1; }

mkdir -p "$SCRATCH_DIR"
echo "Scratch dir : $SCRATCH_DIR"
echo "S3 dest     : $S3_DEST"
echo "Free space  :"
df -h "$SCRATCH_DIR" | sed 's/^/  /'
echo

remote_size() {
  # Content-Length of the remote file, or empty if unavailable
  curl -sSL -I --max-time 60 "$1" \
    | tr -d '\r' \
    | awk 'tolower($1)=="content-length:"{print $2}' \
    | tail -n1
}

s3_size() {
  # Size of the object in S3, or empty if it doesn't exist
  aws s3api head-object --bucket "$1" --key "$2" --query 'ContentLength' --output text 2>/dev/null || true
}

# --- Main loop --------------------------------------------------------------

for entry in "${FILES[@]}"; do
  url_name="${entry%%|*}"
  local_name="${entry##*|}"
  url="$SRC_BASE/$url_name"
  local_path="$SCRATCH_DIR/$local_name"

  # S3 key = everything after s3://bucket/
  bucket="$(echo "$S3_DEST" | sed -E 's#^s3://([^/]+)/.*#\1#')"
  prefix="$(echo "$S3_DEST" | sed -E 's#^s3://[^/]+/(.*)#\1#')"
  key="$prefix/$local_name"

  echo "================================================================"
  echo "File: $local_name"
  echo "  URL: $url"
  echo "  S3 : s3://$bucket/$key"

  expected="$(remote_size "$url")"
  if [[ -z "$expected" ]]; then
    echo "  WARN: could not read remote Content-Length; proceeding without size checks"
  else
    echo "  Remote size: $expected bytes ($(awk -v b="$expected" 'BEGIN{printf "%.1f GB", b/1024/1024/1024}'))"
  fi

  # Skip if already uploaded with the right size
  existing="$(s3_size "$bucket" "$key")"
  if [[ -n "$existing" && "$existing" != "None" ]]; then
    if [[ -z "$expected" || "$existing" == "$expected" ]]; then
      echo "  SKIP: already in S3 (size $existing bytes)"
      continue
    else
      echo "  Present in S3 but size mismatch ($existing != $expected); re-uploading"
    fi
  fi

  # Download with resume + aggressive retry on connection resets
  echo "  Downloading (resumable)..."
  curl -L -C - \
    --retry 100 --retry-all-errors --retry-delay 10 \
    --connect-timeout 60 \
    -o "$local_path" \
    "$url"

  # Verify size before uploading
  if [[ -n "$expected" ]]; then
    actual="$(stat -c%s "$local_path" 2>/dev/null || stat -f%z "$local_path")"
    if [[ "$actual" != "$expected" ]]; then
      echo "  ERROR: downloaded size $actual != expected $expected. Re-run to resume." >&2
      exit 1
    fi
    echo "  Size verified: $actual bytes"
  fi

  # Upload to S3
  echo "  Uploading to s3://$bucket/$key ..."
  aws s3 cp "$local_path" "s3://$bucket/$key"

  # Free scratch space before the next (much larger) file
  echo "  Removing local copy to free space"
  rm -f "$local_path"
  echo "  Done: $local_name"
done

echo "================================================================"
echo "All files mirrored to $S3_DEST"
echo
echo "Next: update the input_raw / input_flat_files URLs in"
echo "  scripts/create_resources/spatial/process_bruker_cosmx_nebius.sh"
echo "to point at the S3 paths, e.g.:"
echo "  input_raw:        $S3_DEST/HalfBrain.zip"
echo "  input_flat_files: $S3_DEST/Half Brain simple files.zip"
echo "  input_raw:        $S3_DEST/NormalLiverFiles.zip"
