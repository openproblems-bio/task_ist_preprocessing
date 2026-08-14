#!/usr/bin/env bash
# fetch_results.sh — copy a Seqera/nebius run's published result files off the
# cluster scratch PVC to a local folder, keyed off the run's WORKFLOW ID.
#
# A benchmark run publishes its outputs to `publish_dir` on the shared scratch PVC
# (tower-scratch, mounted at /scratch), NOT to S3 and NOT anywhere the launch host
# can see. So to get them locally we: (1) ask Seqera for the run's publish_dir, and
# (2) `kubectl cp` the files out of a pod that mounts that PVC.
#
# Usage:
#   fetch_results.sh <workflow_id> [<dest_dir>] [opts]
#
#   <workflow_id>  the Seqera run id (e.g. 3JG1eHsWrAeSse), as recorded in <stage>_runs.csv.
#   <dest_dir>     local folder to copy into (created if absent). Default:
#                  $TXSIM_RESULTS_DIR/runs/<publish-basename> (see --base). Giving an
#                  explicit folder is the normal case (e.g. .../segmentation/cellpose4).
#
#   --all          copy the ENTIRE publish_dir (all files + subdirs, structure preserved),
#                  not just top-level *.yaml (the default). Includes trace.txt etc.
#   --glob <pat>   copy top-level files matching this shell glob instead of '*.yaml'
#                  (e.g. --glob 'score_uns.yaml' for just the scores). Ignored with --all.
#   --workspace <id>  Seqera workspace numeric id (default $TOWER_WORKSPACE_ID or 167877437119966).
#   --pod <name>   pod that mounts the scratch PVC to copy from (default $SCRATCH_POD or 'debug').
#   --namespace <ns>  namespace of that pod / PVC (default $SCRATCH_NS or 'tower-nf').
#   --pvc <name>   PVC to mount if we have to spawn our own pod (default $SCRATCH_PVC or 'tower-scratch').
#   --mount <path> mount path for a spawned pod (default $SCRATCH_MOUNT or '/scratch').
#   --keep-pod     if we spawned an ephemeral pod, leave it running (default: delete it).
#
# If the chosen pod is absent or can't see publish_dir, an ephemeral pod that mounts
# the PVC is spawned (using the process_dataset image, resolved from _viash.yaml), used,
# and torn down (unless --keep-pod). The long-lived `debug` pod is the fast path.
set -uo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || REPO_ROOT="$PWD"
cd "$REPO_ROOT" 2>/dev/null || true

command -v tw      >/dev/null 2>&1 || { echo "ERROR: 'tw' (Seqera CLI) not found on PATH"; exit 2; }
command -v kubectl >/dev/null 2>&1 || { echo "ERROR: 'kubectl' not found on PATH"; exit 2; }

ID="" ; DEST="" ; MODE="yaml" ; GLOB="*.yaml" ; KEEP_POD=0
WS="${TOWER_WORKSPACE_ID:-167877437119966}"
POD="${SCRATCH_POD:-debug}"
NS="${SCRATCH_NS:-tower-nf}"
PVC="${SCRATCH_PVC:-tower-scratch}"
MOUNT="${SCRATCH_MOUNT:-/scratch}"
BASE="${TXSIM_RESULTS_DIR:-$HOME/projects/txsim_results/param_sweep_tests}"

while [ $# -gt 0 ]; do
  case "$1" in
    --all)       MODE="all"; shift;;
    --glob)      GLOB="${2:?}"; shift 2;;
    --workspace) WS="${2:?}"; shift 2;;
    --pod)       POD="${2:?}"; shift 2;;
    --namespace) NS="${2:?}"; shift 2;;
    --pvc)       PVC="${2:?}"; shift 2;;
    --mount)     MOUNT="${2:?}"; shift 2;;
    --base)      BASE="${2:?}"; shift 2;;
    --keep-pod)  KEEP_POD=1; shift;;
    -h|--help)   sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    -*)          echo "unknown option: $1"; exit 2;;
    *)           if [ -z "$ID" ]; then ID="$1"; elif [ -z "$DEST" ]; then DEST="$1";
                 else echo "unexpected arg: $1"; exit 2; fi; shift;;
  esac
done
[ -z "$ID" ] && { echo "usage: fetch_results.sh <workflow_id> [<dest_dir>] [opts]"; exit 2; }
export TOWER_WORKSPACE_ID="$WS"

# ---- resolve publish_dir from Seqera ----
echo "Resolving publish_dir for run $ID (workspace $WS) ..."
params="$(tw runs view -i "$ID" --params 2>&1)" \
  || { echo "$params"; echo "ERROR: could not fetch run $ID (auth? wrong id/workspace?)"; exit 1; }
PUBDIR="$(printf '%s\n' "$params" \
  | sed -nE 's/.*"?publish_dir"?[[:space:]]*[:=][[:space:]]*"?([^",]+)"?.*/\1/p' \
  | head -1 | sed -E 's/[[:space:]]+$//')"
[ -z "$PUBDIR" ] && { echo "ERROR: no publish_dir in the run's params. Params were:"; printf '%s\n' "$params" | head -40; exit 1; }
echo "  publish_dir: $PUBDIR"

# ---- default dest = <base>/runs/<publish-basename> ----
if [ -z "$DEST" ]; then
  DEST="$BASE/runs/$(basename "$PUBDIR")"
  echo "  (no dest given) -> $DEST"
fi

# ---- pick / spawn a pod that mounts the PVC and can see publish_dir ----
SPAWNED=0
pod_sees_dir() { kubectl exec -n "$NS" "$1" -- test -d "$PUBDIR" 2>/dev/null; }
pod_running()  { [ "$(kubectl get pod "$1" -n "$NS" -o jsonpath='{.status.phase}' 2>/dev/null)" = Running ]; }

if pod_running "$POD" && pod_sees_dir "$POD"; then
  echo "Using existing pod '$POD' (namespace $NS)."
else
  if pod_running "$POD"; then
    echo "Pod '$POD' is up but does not see $PUBDIR — spawning an ephemeral pod that mounts $PVC."
  else
    echo "Pod '$POD' not available — spawning an ephemeral pod that mounts $PVC."
  fi
  POD="fetch-results-${ID}"
  # image ref from _viash.yaml (same resolution the cluster pulls with); guaranteed pullable.
  ident() { grep -E "$1" _viash.yaml 2>/dev/null | head -1 | sed -E 's/.*:[[:space:]]*//; s/[[:space:]]*$//'; }
  ORG="$(ident '^organization:')"; PROJ="$(ident '^name:')"
  REG="$(ident 'docker_registry:')"; REG="${REG:-ghcr.io}"
  IMAGE="${FETCH_IMAGE:-$REG/$ORG/$PROJ/data_processors/process_dataset:build_main}"
  echo "  image: $IMAGE"
  cat <<YAML | kubectl apply -f - 2>&1 | tail -1
apiVersion: v1
kind: Pod
metadata:
  name: $POD
  namespace: $NS
spec:
  restartPolicy: Never
  containers:
  - name: main
    image: $IMAGE
    command: ["sleep", "infinity"]
    volumeMounts:
    - name: scratch
      mountPath: $MOUNT
  volumes:
  - name: scratch
    persistentVolumeClaim:
      claimName: $PVC
YAML
  SPAWNED=1
  echo "  waiting for Ready (image may be multi-GB; cached nodes are instant)..."
  kubectl wait --for=condition=Ready "pod/$POD" -n "$NS" --timeout=600s 2>&1 | tail -1
  pod_sees_dir "$POD" || { echo "ERROR: spawned pod still can't see $PUBDIR (wrong PVC/mount?)."; \
    [ "$KEEP_POD" -eq 0 ] && kubectl delete pod "$POD" -n "$NS" --wait=false >/dev/null 2>&1; exit 1; }
fi

cleanup() { [ "$SPAWNED" -eq 1 ] && [ "$KEEP_POD" -eq 0 ] && \
  { echo "tearing down ephemeral pod '$POD' ..."; kubectl delete pod "$POD" -n "$NS" --wait=false >/dev/null 2>&1; }; }
trap cleanup EXIT

mkdir -p "$DEST"

# ---- copy ----
if [ "$MODE" = all ]; then
  echo "Copying ENTIRE publish_dir -> $DEST ..."
  # kubectl cp <pod>:<dir> <dest> lands the dir's CONTENTS directly in <dest> (which we
  # already mkdir'd), preserving any subdirs — not nested under <dest>/<basename>.
  kubectl cp "$NS/$POD:$PUBDIR" "$DEST" 2>&1 | grep -v 'Removing leading' || true
else
  files="$(kubectl exec -n "$NS" "$POD" -- bash -c "cd '$PUBDIR' 2>/dev/null && ls -1 $GLOB 2>/dev/null")"
  [ -z "$files" ] && { echo "ERROR: no files matching '$GLOB' in $PUBDIR"; \
    echo "  (dir listing:)"; kubectl exec -n "$NS" "$POD" -- bash -c "ls -1 '$PUBDIR'" 2>&1 | sed 's/^/    /'; exit 1; }
  echo "Copying files matching '$GLOB' -> $DEST ..."
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    echo "  $f"
    kubectl cp "$NS/$POD:$PUBDIR/$f" "$DEST/$f" 2>&1 | grep -v 'Removing leading' || true
  done <<< "$files"
fi

# ---- report ----
echo
echo "Done -> $DEST"
ls -la "$DEST"
