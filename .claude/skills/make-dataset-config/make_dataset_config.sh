#!/usr/bin/env bash
# make_dataset_config.sh — create a Dataset_report dataset-paths YAML for a combined
# dataset that lives on the cluster scratch PVC, and DOWNLOAD its raw single-cell
# reference into raw/ if that .h5ad is missing on the pod. Writes results/config/<id>.yaml
# with an id that matches the dataset's id on the pod exactly.
#
# The raw single-cell reference is the combine step's input_sc. For standard datasets
# that path is `$input_dir/<subpath>` where `$input_dir=s3://openproblems-data/resources/
# datasets`, so it mirrors 1:1 onto the PVC (`raw/<subpath>`) and onto the public HTTPS
# endpoint (https://openproblems-data.s3.amazonaws.com/resources/datasets/<subpath>). The
# helper resolves that mapping straight from the combine scripts; the spatial + panel-subset
# SC come from the combined outputs (output_sp.zarr / output_sc.h5ad) already on the PVC.
#
# Usage:
#   make_dataset_config.sh <dataset_id> [--sc-source <subpath>] [--sc-raw <raw-rel-path>] [--out <file>]
#
#   <dataset_id>        the combined dataset's id = its dir under datasets/ on the pod
#                       (e.g. 2024_10x_human_liver_xenium_combined, or ".../rep1").
#   --sc-source <sub>   override the raw-SC S3 subpath under resources/datasets/ (download-if-missing).
#                       Use when the combine-script auto-resolution can't (e.g. a $sc_ref var).
#   --sc-raw <path>     the raw SC is ALREADY on the PVC at this raw-relative path (a custom
#                       reference, e.g. raw/ltx_human_lung_sc/.../dataset.h5ad); no download.
#   --out <file>        output path (default results/config/<id-with-slashes->_>.yaml).
#
# Env overrides (defaults target this repo's Nebius setup):
#   SCRATCH_POD (debug)  SCRATCH_NS (tower-nf)  SCRATCH_MOUNT (/scratch)
#   DATASET_ROOT ($SCRATCH_MOUNT/task_ist_preprocessing)  CONFIG_DIR (results/config)
#   S3_BASE (https://openproblems-data.s3.amazonaws.com/resources/datasets)
set -uo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "ERROR: not in a git repo"; exit 2; }
cd "$REPO_ROOT" || exit 2

ID=""; SC_SRC=""; SC_RAW=""; OUT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --sc-source) SC_SRC="${2:?}"; shift 2;;
    --sc-raw)    SC_RAW="${2:?}"; shift 2;;
    --out)       OUT="${2:?}";    shift 2;;
    -h|--help)   grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 0;;
    -*)          echo "unknown flag: $1"; exit 2;;
    *)           [ -z "$ID" ] && ID="$1"; shift;;
  esac
done
[ -z "$ID" ] && { grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 2; }

POD="${SCRATCH_POD:-debug}"
NS="${SCRATCH_NS:-tower-nf}"
MOUNT="${SCRATCH_MOUNT:-/scratch}"
ROOT="${DATASET_ROOT:-$MOUNT/task_ist_preprocessing}"
S3_BASE="${S3_BASE:-https://openproblems-data.s3.amazonaws.com/resources/datasets}"
CONFIG_DIR="${CONFIG_DIR:-results/config}"

xp() { kubectl exec -n "$NS" "$POD" -- sh -c "$1"; }

kubectl get pod "$POD" -n "$NS" >/dev/null 2>&1 || {
  echo "ERROR: pod '$POD' not found in namespace '$NS'. Point SCRATCH_POD at any pod that"
  echo "       mounts the scratch PVC (with curl + python3/h5py), or bring one up via the"
  echo "       render-qmd-pod skill (\$S up) and set SCRATCH_POD=render-qmd."
  exit 2
}

echo "Dataset:  $ID"
echo "Pod:      $POD (ns $NS)   root: $ROOT"

# ── 1. the combined dataset dir must exist with its outputs ──────────────────
DDIR="datasets/$ID"
have="$(xp "cd '$ROOT/$DDIR' 2>/dev/null && ls -d output_sp.zarr output_sc.h5ad 2>/dev/null | tr '\n' ' '")"
case "$have" in
  *output_sp.zarr*output_sc.h5ad*|*output_sc.h5ad*output_sp.zarr*) : ;;
  *) echo "ERROR: $ROOT/$DDIR is missing output_sp.zarr and/or output_sc.h5ad (got: '${have:-none}')"
     echo "       Is the dataset id correct? List: kubectl exec -n $NS $POD -- ls $ROOT/datasets"
     exit 3;;
esac
echo "Combined: $DDIR  (output_sp.zarr + output_sc.h5ad present)"

# ── 2. resolve the raw single-cell reference: RAW_REL (+ optional URL) ────────
RAW_REL=""; URL=""
if [ -n "$SC_RAW" ]; then
  RAW_REL="$SC_RAW"                                   # already-on-PVC custom ref, no download
elif [ -n "$SC_SRC" ]; then
  RAW_REL="raw/$SC_SRC"; URL="$S3_BASE/$SC_SRC"       # explicit S3 subpath
else
  # auto: pull input_sc for this id from the combine scripts (first match wins).
  raw_sc_line="$(cat scripts/create_resources/combine/*.sh 2>/dev/null | awk -v id="$ID" '
    index($0, "id: \"" id "\"") {f=1}
    f && /input_sc:/ {print; exit}')"
  val="$(printf '%s' "$raw_sc_line" | sed -E 's/.*input_sc:[[:space:]]*"([^"]*)".*/\1/')"
  case "$val" in
    '$input_dir/'*|'${input_dir}/'*)
      sub="${val#*input_dir\}/}"; sub="${sub#*input_dir/}"
      RAW_REL="raw/$sub"; URL="$S3_BASE/$sub";;
    *)
      echo "ERROR: could not auto-resolve the raw SC reference for '$ID'."
      echo "       input_sc found: '${val:-<none>}'"
      echo "       Pass --sc-source <subpath-under-resources/datasets> (to download from S3),"
      echo "       or --sc-raw <raw-relative-path-already-on-the-PVC>."
      exit 3;;
  esac
fi
echo "Raw SC:   $RAW_REL"
[ -n "$URL" ] && echo "          (S3: $URL)"

# ── 3. download the raw SC .h5ad to raw/ if missing (idempotent) ─────────────
if xp "test -f '$ROOT/$RAW_REL'"; then
  echo "          already present ($(xp "du -h '$ROOT/$RAW_REL' 2>/dev/null | cut -f1"))"
elif [ -n "$URL" ]; then
  echo "          missing -> downloading from S3 ..."
  xp "set -e; d=\$(dirname '$ROOT/$RAW_REL'); mkdir -p \"\$d\";
      curl -sSf '$URL' -o '$ROOT/$RAW_REL.part' && mv '$ROOT/$RAW_REL.part' '$ROOT/$RAW_REL'" \
    || { echo "ERROR: download failed from $URL"; exit 4; }
  echo "          downloaded ($(xp "du -h '$ROOT/$RAW_REL' 2>/dev/null | cut -f1"))"
else
  echo "ERROR: raw SC ref missing at $ROOT/$RAW_REL and no download URL (custom --sc-raw)."
  exit 4
fi

# ── 4. best-effort sanity check that the ref has what the report needs ───────
kubectl exec -i -n "$NS" "$POD" -- python - "$ROOT/$RAW_REL" <<'PY' 2>/dev/null || echo "          (skipped structure check: no python/h5py in $POD)"
import sys, h5py
p = sys.argv[1]
with h5py.File(p, "r") as f:
    obsm   = set(f["obsm"].keys())   if "obsm"   in f else set()
    layers = set(f["layers"].keys()) if "layers" in f else set()
    obs    = set(f["obs"].keys())    if "obs"    in f else set()
    miss = [n for n,ok in [("obsm/X_pca","X_pca" in obsm),
                            ("layers/normalized","normalized" in layers),
                            ("obs/cell_type", any(k=="cell_type" for k in obs))] if not ok]
    print("          structure:", "OK (X_pca, normalized, cell_type)" if not miss
          else "WARNING missing " + ", ".join(miss))
PY

# ── 5. write the config YAML ─────────────────────────────────────────────────
SLUG="$(printf '%s' "$ID" | tr '/' '_')"
OUT="${OUT:-$CONFIG_DIR/$SLUG.yaml}"
mkdir -p "$(dirname "$OUT")"
SRC_COMMENT="already staged on the PVC (custom reference)"; [ -n "$URL" ] && SRC_COMMENT="$URL"
cat > "$OUT" <<YAML
# Input paths for Dataset_report.qmd — Nebius k8s (tower-scratch PVC).
# Generated by .claude/skills/make-dataset-config/make_dataset_config.sh
#
# Dataset id — matches the dataset dir / dataset_id on the pod exactly:
id: $ID
#
# Paths are resolved relative to the resources root:
#   \$RESOURCES_DIR, else the first existing entry in resources_roots.
# Override the location of this config file with \$DATASET_CONFIG.

resources_roots:
  - $ROOT                                                # Nebius k8s: tower-scratch PVC mount
  - ~/projects/task_ist_preprocessing/resources_test     # local checkout fallback

# Full single-cell reference (raw source): processed (counts + normalized layers, X_pca),
# full gene set. This is the combine step's input_sc; downloaded to raw/ if missing from:
#   $SRC_COMMENT
scrnaseq_raw: $RAW_REL

# Combined process_datasets output for this dataset (spatial + panel-subset SC).
dataset:
  dir: $DDIR
  scrnaseq_panel: output_sc.h5ad   # reference subsetted to the spatial panel
  spatial: output_sp.zarr          # combined IST input (image, cell_labels, transcripts)
YAML
echo "Wrote:    $OUT"
