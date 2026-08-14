#!/bin/bash

# Upload a 2D-only copy of the Kuppe kidney MERFISH (Vizgen MERSCOPE) raw data to S3.
#
# The raw region_0 folders carry 7 focal z-planes per channel (mosaic_*_z0..z6.tif),
# ~500-820 GB per sample. The repo's vizgen_merscope loader reads ONLY z-plane 3
# (spatialdata_io.merscope(..., z_layers=3), which globs the images dir for stains and
# opens mosaic_<stain>_z3.tif). So we mirror a z3-only copy: lossless for this 2D
# pipeline, ~79 GB (cancer) + ~135 GB (healthy).
#
# Kept per sample (region_0/): mosaic_<5 stains>_z3.tif, images/manifest.json,
# images/micron_to_mosaic_pixel_transform.csv, boundary_mask_images/boundaries_z3.tif,
# detected_transcripts.csv (all molecules), cell_by_gene.csv, cell_metadata.csv,
# cell_boundaries.parquet, summary.png.
# Dropped: every mosaic_*_z{0,1,2,4,5,6}.tif and boundaries_z{0,1,2,4,5,6}.tif.
#
# Idempotent + resumable: each file is skipped if already in S3 with a matching size,
# so just re-run if the connection drops. aws s3 cp handles multipart for the big tifs.
#
# Usage:
#   SRC_BASE=/Volumes/SeagateHHD/Kuppe_kidney_merfish/240824_download/merfish \
#   S3_DEST=s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney/merfish \
#   AWS_PROFILE=op ./upload_kuppe_kidney_merfish_2d.sh

set -euo pipefail

SRC_BASE="${SRC_BASE:-/Volumes/SeagateHHD/Kuppe_kidney_merfish/240824_download/merfish}"
S3_DEST="${S3_DEST:-s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney/merfish}"
AWS_PROFILE="${AWS_PROFILE:-op}"
SAMPLES=(${SAMPLES:-cancer healthy})
Z="${Z:-3}"   # which single z-plane to keep

export AWS_PROFILE
command -v aws >/dev/null || { echo "ERROR: aws CLI not found" >&2; exit 1; }

bucket="$(echo "$S3_DEST" | sed -E 's#^s3://([^/]+)/.*#\1#')"
prefix="$(echo "$S3_DEST" | sed -E 's#^s3://[^/]+/(.*)#\1#')"

file_size() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }
s3_size()   { aws s3api head-object --bucket "$1" --key "$2" --query 'ContentLength' --output text 2>/dev/null || true; }

# Build the wanted-file list for a sample's region_0 (relative paths under region_0/).
wanted_rel() {
  local r="$SRC_BASE/$1/region_0"
  ( cd "$r" && \
    ls images/mosaic_*_z${Z}.tif \
       images/manifest.json \
       images/micron_to_mosaic_pixel_transform.csv \
       boundary_mask_images/boundaries_z${Z}.tif \
       detected_transcripts.csv cell_by_gene.csv cell_metadata.csv \
       cell_boundaries.parquet summary.png 2>/dev/null )
}

total_files=0; uploaded=0; skipped=0
for sample in "${SAMPLES[@]}"; do
  echo "================================================================"
  echo "Sample: $sample   ($SRC_BASE/$sample/region_0)"
  while IFS= read -r rel; do
    [ -z "$rel" ] && continue
    src="$SRC_BASE/$sample/region_0/$rel"
    key="$prefix/$sample/region_0/$rel"
    total_files=$((total_files+1))
    local_sz="$(file_size "$src")"
    remote_sz="$(s3_size "$bucket" "$key")"
    if [[ -n "$remote_sz" && "$remote_sz" != "None" && "$remote_sz" == "$local_sz" ]]; then
      printf '  SKIP  %-55s (%s bytes, already in S3)\n' "$rel" "$local_sz"
      skipped=$((skipped+1))
      continue
    fi
    printf '  UP    %-55s (%.2f GB)\n' "$rel" "$(awk -v b="$local_sz" 'BEGIN{print b/1e9}')"
    aws s3 cp "$src" "s3://$bucket/$key" --only-show-errors
    uploaded=$((uploaded+1))
  done < <(wanted_rel "$sample")
done

echo "================================================================"
echo "Done. files=$total_files uploaded=$uploaded skipped=$skipped"
echo "Loader --input paths:"
for sample in "${SAMPLES[@]}"; do
  echo "  $S3_DEST/$sample/region_0"
done
