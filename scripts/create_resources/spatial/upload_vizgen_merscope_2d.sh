#!/bin/bash

# Mirror a z3-only (2D) copy of the Vizgen MERSCOPE FFPE showcase raw data from Google
# Cloud Storage to S3, so the Nebius/Seqera benchmark can stage its inputs without any
# Google Cloud credentials.
#
# WHY THIS EXISTS
#   process_vizgen_merscope_nebius.sh used to point `input:` directly at
#   gs://vz-ffpe-showcase/<Patient>. That bucket is Vizgen's access-controlled
#   data-release bucket (NOT public: an unauthenticated GET returns
#   "Anonymous caller does not have storage.objects.get access ..."). The Nebius compute
#   env has no GCP credentials, so Nextflow staged the input as an anonymous caller and
#   the run died before the loader ran. Mirroring the data to s3://openproblems-data
#   (which the Nebius env already reads for every other spatial dataset) removes the
#   GCP-auth dependency from every future run.
#
# WHAT IS KEPT (lossless for this 2D pipeline)
#   The vizgen_merscope loader calls spatialdata_io.merscope(..., z_layers=3), i.e. it
#   only ever reads z-plane 3 of the mosaic images. So we mirror ONLY *_z3.tif and drop
#   every *_z{0,1,2,4,5,6}.tif (7 focal planes -> 1). Everything else is copied verbatim:
#     - cell_boundaries/*.hdf5      (OLD showcase format; the loader rebuilds
#                                    cell_boundaries.parquet from these at run time)
#     - images/mosaic_DAPI_z3.tif   (DAPI only, by default — the loader keeps ONLY the DAPI
#                                    channel: `sdata["morphology_mip"].sel(c=["DAPI"])`, so
#                                    dropping PolyT/Cellbound stains is lossless for this
#                                    pipeline and cuts the image bytes to ~1/5. If the
#                                    keep-more-stains TODO in the loader is ever done, re-run
#                                    with KEEP_ONLY_DAPI=0 to mirror all stains at z3.)
#     - images/micron_to_mosaic_pixel_transform.csv, images/manifest.json
#     - cell_by_gene.csv, cell_metadata.csv, detected_transcripts.csv
#   KEEP_ONLY_DAPI defaults to 1 (DAPI-only). Set KEEP_ONLY_DAPI=0 to keep all stains at z3.
#
# PARALLEL
#   Each patient ships ~1100-2500 tiny cell_boundaries/*.hdf5 files (~15k total across the 7
#   patients). Streaming them one-at-a-time is dominated by per-file spawn/handshake overhead
#   (~15 s/file => ~2-3 DAYS). This script runs the transfers through a pool of $JOBS parallel
#   workers (default 16), which hides that overhead and cuts the whole mirror to a few hours.
#   Big z3 tifs / detected_transcripts.csv are bandwidth-bound, so a handful running
#   concurrently just keeps the uplink saturated. Tune with JOBS=<n> (e.g. 32 for the
#   mostly-small-file tail).
#
# REQUIREMENTS (run this on a machine that is authenticated to BOTH clouds)
#   - gcloud CLI, authenticated with an account granted access to gs://vz-ffpe-showcase
#     (Vizgen data-release-program access — request it via https://info.vizgen.com/ffpe-showcase).
#   - aws CLI, with a profile that can write s3://openproblems-data (default: profile "op").
#   The transfer STREAMS each object gcloud->aws (no local disk staging).
#
# IDEMPOTENT + RESUMABLE
#   Objects already in S3 with a matching byte size are skipped (one bulk `aws s3 ls` up
#   front, so resuming does not cost a HEAD per object). A partially-copied file is aborted
#   by aws and reads back as absent, so it is simply re-sent on the next run. Just re-run.
#
# Usage:
#   AWS_PROFILE=op ./upload_vizgen_merscope_2d.sh
#   # more workers for the small-file tail:
#   JOBS=32 AWS_PROFILE=op ./upload_vizgen_merscope_2d.sh
#   # or override anything:
#   GCS_BASE=gs://vz-ffpe-showcase \
#   S3_DEST=s3://openproblems-data/resources/raw_data/txSim_custom/vizgen_merscope \
#   SAMPLES="HumanBreastCancerPatient1 HumanLungCancerPatient1" \
#   JOBS=24 AWS_PROFILE=op ./upload_vizgen_merscope_2d.sh

set -euo pipefail

GCS_BASE="${GCS_BASE:-gs://vz-ffpe-showcase}"
S3_DEST="${S3_DEST:-s3://openproblems-data/resources/raw_data/txSim_custom/vizgen_merscope}"
AWS_PROFILE="${AWS_PROFILE:-op}"
Z="${Z:-3}"                             # which single z-plane to keep
KEEP_ONLY_DAPI="${KEEP_ONLY_DAPI:-1}"   # 1 = DAPI-only (default, lossless for this loader); 0 = all stains at z3
JOBS="${JOBS:-12}"                      # number of concurrent gcloud->aws transfers (each = 2 TLS
                                        # conns; >~16 tends to trigger transient TLS/read-timeout
                                        # errors, which the worker retries with backoff anyway)

# The 7 patients currently active (uncommented) in process_vizgen_merscope_nebius.sh.
# Additional patients (melanoma/ovarian/prostate/uterine) are commented-out there and in
# the s3 sibling process_vizgen_merscope.sh — add their folder names here to mirror them.
SAMPLES=(${SAMPLES:-\
  HumanBreastCancerPatient1 \
  HumanLiverCancerPatient1 \
  HumanLiverCancerPatient2 \
  HumanLungCancerPatient1 \
  HumanLungCancerPatient2 \
  HumanColonCancerPatient1 \
  HumanColonCancerPatient2})

export AWS_PROFILE
command -v gcloud >/dev/null || { echo "ERROR: gcloud CLI not found (needed to read gs://vz-ffpe-showcase)" >&2; exit 1; }
command -v aws    >/dev/null || { echo "ERROR: aws CLI not found (needed to write $S3_DEST)" >&2; exit 1; }

gcs_bucket="$(echo "$GCS_BASE" | sed -E 's#^gs://([^/]+).*#\1#')"
s3_bucket="$(echo "$S3_DEST"  | sed -E 's#^s3://([^/]+)/.*#\1#')"
s3_prefix="$(echo "$S3_DEST"  | sed -E 's#^s3://[^/]+/(.*)#\1#')"
export gcs_bucket s3_bucket s3_prefix

# Should this object be dropped from the mirror?
drop_object() {
  local rel="$1" base="${1##*/}"
  # keep only z-plane $Z: drop any *_z<other>.tif (matches mosaic_*_z*.tif and boundaries_z*.tif)
  if [[ "$base" =~ _z([0-9]+)\.tif$ && "${BASH_REMATCH[1]}" != "$Z" ]]; then return 0; fi
  # optional: drop non-DAPI stains entirely
  if [[ "$KEEP_ONLY_DAPI" == "1" && "$base" =~ ^mosaic_ && ! "$base" =~ ^mosaic_DAPI_ ]]; then return 0; fi
  return 1
}

# Worker: stream one "<size>|<rel>" record gcloud->aws. Exported so `xargs -P` can run it in
# parallel `bash -c` subshells. Kept free of bash-4 features so it works on the stock macOS
# /bin/bash 3.2. The already-uploaded files are filtered out up front (see below), so the
# worker just transfers. A failed transfer is logged, not fatal: the batch keeps going and the
# file (absent, since aws aborts a partial multipart) is retried on the next resumable run.
transfer_one() {
  local rec="$1"
  local size="${rec%%|*}"
  local rel="${rec#*|}"
  local attempt=0 max=4
  # Retry transient GCS/S3 errors (read timeouts, SSL/TLS handshake failures) that show up under
  # high concurrency. pipefail (set by the caller) makes a gcloud-side failure fail the pipe too.
  while :; do
    attempt=$((attempt+1))
    if gcloud storage cat "gs://$gcs_bucket/$rel" 2>/dev/null \
         | aws s3 cp - "s3://$s3_bucket/$s3_prefix/$rel" --expected-size "$size" --only-show-errors 2>/dev/null; then
      printf '  OK    %s (%s B)\n' "$rel" "$size"
      return 0
    fi
    if [ "$attempt" -ge "$max" ]; then
      printf '  FAIL  %s (after %s attempts; retry on re-run)\n' "$rel" "$attempt" >&2
      return 0
    fi
    sleep $((attempt * 3))   # linear backoff: 3s, 6s, 9s
  done
}
export -f transfer_one

present="$(mktemp -t viz_present)"
wanted_raw="$(mktemp -t viz_wanted_raw)"
wanted="$(mktemp -t viz_wanted)"
worklist="$(mktemp -t viz_worklist)"
trap 'rm -f "$present" "$wanted_raw" "$wanted" "$worklist"' EXIT

# 1) One bulk listing of what is already in S3 -> "<rel>\t<size>" (sorted). One API call for the
#    whole resume state, instead of a HEAD per object.
echo "Listing existing objects under $S3_DEST ..."
aws s3 ls "s3://$s3_bucket/$s3_prefix/" --recursive 2>/dev/null \
  | awk -v p="$s3_prefix/" 'BEGIN{OFS="\t"} {k=$4; sub("^"p,"",k); if(k!="") print k,$3}' \
  | sort > "$present"
echo "  $(wc -l < "$present" | tr -d ' ') objects already in S3."

# 2) Enumerate GCS, apply the keep-rules -> "<rel>\t<size>". Written to a file (not piped) so
#    the counters stay in this shell and the progress echoes go to stderr, not into the data.
total=0; dropped=0
for sample in "${SAMPLES[@]}"; do
  echo "Enumerating $GCS_BASE/$sample ..." >&2
  while IFS= read -r line; do
    [[ "$line" == *"gs://"* ]] || continue           # skip the TOTAL: summary line
    url="$(awk '{print $NF}' <<<"$line")"
    size="$(awk '{print $1}'  <<<"$line")"
    [[ "$url" == */ ]] && continue                   # skip directory placeholders
    [[ "$size" =~ ^[0-9]+$ ]] || continue            # skip anything without a numeric size
    rel="${url#gs://$gcs_bucket/}"
    total=$((total+1))
    if drop_object "$rel"; then dropped=$((dropped+1)); continue; fi
    printf '%s\t%s\n' "$rel" "$size" >> "$wanted_raw"
  done < <(gcloud storage ls -l "$GCS_BASE/$sample/**")
done
sort "$wanted_raw" > "$wanted"

# 3) Work list = wanted objects not already in S3 at the same size (comm on the sorted files).
comm -23 "$wanted" "$present" | awk -F'\t' '{print $2"|"$1}' > "$worklist"
queued="$(wc -l < "$worklist" | tr -d ' ')"
echo "================================================================"
echo "Enumerated $total objects: $dropped dropped, $(( $(wc -l < "$wanted" | tr -d ' ') - queued )) already in S3, $queued to transfer."
echo "Transferring with $JOBS parallel workers ..."

# 4) Run the transfers in parallel. -I REC feeds one record per worker; -P runs $JOBS at once.
#    Use the SAME bash that is running this script ($BASH) so the exported transfer_one
#    imports correctly regardless of which bash is first on PATH.
if [[ "$queued" -gt 0 ]]; then
  xargs -P "$JOBS" -I REC "${BASH:-bash}" -c 'set -o pipefail; transfer_one "$@"' _ REC < "$worklist"
fi

echo "================================================================"
echo "Done. objects=$total dropped=$dropped queued=$queued (JOBS=$JOBS)"
echo "Loader --input paths (use these in process_vizgen_merscope_nebius.sh):"
for sample in "${SAMPLES[@]}"; do
  echo "  $S3_DEST/$sample"
done
