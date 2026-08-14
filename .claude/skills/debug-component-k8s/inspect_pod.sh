#!/usr/bin/env bash
# inspect_pod.sh — spin up a throwaway pod running a Viash component's container
# image on the CURRENT kubectl cluster, to troubleshoot its runtime environment
# (imports, dependency versions, native libraries, CLI) WITHOUT a full Nextflow /
# Seqera benchmark run. This is the interactive counterpart to `check-component`
# (which only checks deployment freshness).
#
# The image is resolved from _viash.yaml the same way the cluster pulls it:
#   <docker_registry>/<organization>/<name>/<namespace>/<component>:<tag>
#
# Usage:
#   inspect_pod.sh up    <component> [tag] [--gpu] [--scratch]  # launch pod + wait Ready
#   inspect_pod.sh run   <component> '<bash>'        # run a bash snippet in the pod
#   inspect_pod.sh py    <component> '<python>'      # run python (RAPIDS no-init env preset)
#   inspect_pod.sh down  <component>                 # delete the pod
#   inspect_pod.sh image <component> [tag]           # just print the resolved image ref
#
#   tag defaults to build_main. Add --gpu to schedule on a GPU node (for tests that
#   actually need a device — rmm/cuInit, a real run). Most dep/import/native-lib
#   troubleshooting needs NO gpu, so omit it (faster, schedules anywhere).
#
#   Add --scratch to mount the shared scratch PVC at $SCRATCH_MOUNT so you can inspect a
#   Nextflow work dir's staged inputs/outputs (h5ad/zarr) using the COMPONENT'S OWN image —
#   its real libraries (anndata/spatialdata/...) read those files natively, with no ad-hoc
#   pip installs in a generic pod. PVCs are namespace-scoped, so the pod then runs in the
#   PVC's namespace. You must know the exact work-dir path (an opaque hash) — ASK THE USER
#   for it; it can't be guessed. Cluster specifics default to this repo's Seqera setup and
#   can be overridden via env:
#     SCRATCH_PVC (tower-scratch)   SCRATCH_NS (tower-nf)   SCRATCH_MOUNT (/scratch)
set -uo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "ERROR: not in a git repo"; exit 2; }
cd "$REPO_ROOT" || exit 2

CMD="${1:-}"; ARG="${2:-}"
[ -z "$CMD" ] || [ -z "$ARG" ] && { grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 2; }

# ---- resolve component dir -> namespace/component ----
if [[ -d "src/$ARG" && -f "src/$ARG/config.vsh.yaml" ]]; then
  COMP_DIR="src/$ARG"
else
  HITS="$(find src -type d -name "$ARG" 2>/dev/null | while read -r d; do [ -f "$d/config.vsh.yaml" ] && echo "$d"; done)"
  N="$(printf '%s\n' "$HITS" | grep -c . || true)"
  [ "$N" -eq 0 ] && { echo "ERROR: no component '$ARG' under src/"; exit 2; }
  [ "$N" -gt 1 ] && { echo "ERROR: '$ARG' is ambiguous:"; printf '%s\n' "$HITS"; exit 2; }
  COMP_DIR="$HITS"
fi
COMP="$(basename "$COMP_DIR")"
NS="${COMP_DIR#src/}"; NS="${NS%/$COMP}"

# ---- project identity from _viash.yaml ----
ident() { grep -E "$1" _viash.yaml | head -1 | sed -E 's/.*:[[:space:]]*//; s/[[:space:]]*$//'; }
ORG="$(ident '^organization:')"; PROJ="$(ident '^name:')"
REG="$(ident 'docker_registry:')"; REG="${REG:-ghcr.io}"

# ---- tag / flags ----
TAG="build_main"; GPU=0; SCR=0
for a in "${@:3}"; do
  case "$a" in
    --gpu) GPU=1 ;;
    --scratch) SCR=1 ;;
    -*) : ;;
    *) TAG="$a" ;;
  esac
done
IMG="$REG/$ORG/$PROJ/$NS/$COMP:$TAG"
POD="inspect-$(printf '%s' "$COMP" | tr '_' '-' | tr '[:upper:]' '[:lower:]')"

# ---- scratch PVC config (namespace-scoped; override via env for other clusters) ----
SCRATCH_PVC="${SCRATCH_PVC:-tower-scratch}"
SCRATCH_NS="${SCRATCH_NS:-tower-nf}"
SCRATCH_MOUNT="${SCRATCH_MOUNT:-/scratch}"

# ---- GPU scheduling (override via env for other clusters) ----
# GPU nodes on this Nebius cluster sit in a dedicated node group (see the `gpu` label in
# src/base/labels_nebius.config) and may carry a taint, so a bare `nvidia.com/gpu: 1` limit
# is not enough to place the pod — we also pin the nodeSelector and tolerate any taint.
GPU_NODESELECTOR_KEY="${GPU_NODESELECTOR_KEY:-nebius.com/node-group-id}"
GPU_NODESELECTOR_VAL="${GPU_NODESELECTOR_VAL:-mk8snodegroup-e00t775jb99svb7k5r}"

# Resolve the pod's namespace by name: a plain pod lands in the current ns, but a
# --scratch pod lives in SCRATCH_NS. run/py/down use this so they hit it either way.
# (Guarded array expansion below keeps this working under bash 3.2 + `set -u`.)
pod_ns() { kubectl get pod -A --field-selector "metadata.name=$POD" -o jsonpath='{.items[0].metadata.namespace}' 2>/dev/null; }

case "$CMD" in
  image) echo "$IMG" ;;

  up)
    echo "Component: $COMP ($NS)   Pod: $POD"
    echo "Image:     $IMG"
    UPNS=()
    if [ "$SCR" -eq 1 ]; then
      # PVCs are namespace-scoped, so the pod must live in the PVC's namespace.
      echo "Scratch:   PVC $SCRATCH_PVC -> $SCRATCH_MOUNT   (namespace $SCRATCH_NS)"
      UPNS=(-n "$SCRATCH_NS")
      GPU_BLOCK=""
      [ "$GPU" -eq 1 ] && GPU_BLOCK="    resources:
      limits:
        nvidia.com/gpu: 1"
      # kubectl run can't attach a PVC volume; apply a manifest instead. Same
      # `sleep infinity` rationale as below (don't let the pod complete).
      cat <<YAML | kubectl apply -f - 2>&1 | tail -1
apiVersion: v1
kind: Pod
metadata:
  name: $POD
  namespace: $SCRATCH_NS
spec:
  restartPolicy: Never
  containers:
  - name: main
    image: $IMG
    command: ["sleep", "infinity"]
$GPU_BLOCK
    volumeMounts:
    - name: scratch
      mountPath: $SCRATCH_MOUNT
  volumes:
  - name: scratch
    persistentVolumeClaim:
      claimName: $SCRATCH_PVC
YAML
    elif [ "$GPU" -eq 1 ]; then
      # Request a GPU + pin to the GPU node group and tolerate its taint. NOTE: `kubectl run
      # --limits` was REMOVED in modern kubectl, so build the pod from a manifest (same
      # `sleep infinity` rationale as below). Override GPU_NODESELECTOR_KEY/VAL for other
      # clusters; clear them to rely on the resource limit alone.
      NODESEL=""
      [ -n "$GPU_NODESELECTOR_KEY" ] && NODESEL="  nodeSelector:
    $GPU_NODESELECTOR_KEY: $GPU_NODESELECTOR_VAL"
      cat <<YAML | kubectl apply -f - 2>&1 | tail -1
apiVersion: v1
kind: Pod
metadata:
  name: $POD
spec:
  restartPolicy: Never
$NODESEL
  tolerations:
  - operator: "Exists"
  containers:
  - name: main
    image: $IMG
    command: ["sleep", "infinity"]
    resources:
      limits:
        nvidia.com/gpu: 1
YAML
    else
      # NOTE: use `sleep infinity` — a finite `sleep N` (or the image entrypoint)
      # can let the pod complete (phase Succeeded), and you can't exec a completed pod.
      kubectl run "$POD" --image="$IMG" --restart=Never \
        --command -- sleep infinity 2>&1 | tail -1
    fi
    echo "waiting for Ready (image may be multi-GB; already-cached nodes are instant)..."
    kubectl ${UPNS[@]+"${UPNS[@]}"} wait --for=condition=Ready "pod/$POD" --timeout=600s 2>&1 | tail -1
    kubectl ${UPNS[@]+"${UPNS[@]}"} get pod "$POD" -o wide 2>&1 | tail -1
    ;;

  run)
    SNIP="${3:?need a bash snippet as the 3rd arg}"
    PNS="$(pod_ns)"; NSARG=(); [ -n "$PNS" ] && NSARG=(-n "$PNS")
    kubectl ${NSARG[@]+"${NSARG[@]}"} exec "$POD" -- bash -c "$SNIP" 2>&1
    ;;

  py)
    SNIP="${3:?need a python snippet as the 3rd arg}"
    PNS="$(pod_ns)"; NSARG=(); [ -n "$PNS" ] && NSARG=(-n "$PNS")
    # preset RAPIDS no-init so `import cudf`/segger deps don't abort on a GPU-less pod
    kubectl ${NSARG[@]+"${NSARG[@]}"} exec "$POD" -- bash -c '
      export RAPIDS_NO_INITIALIZE=1 CUDF_NO_INITIALIZE=1 RMM_NO_INITIALIZE=1
      python - <<'"'"'PYEOF'"'"'
'"$SNIP"'
PYEOF' 2>&1
    ;;

  down)
    PNS="$(pod_ns)"; NSARG=(); [ -n "$PNS" ] && NSARG=(-n "$PNS")
    kubectl ${NSARG[@]+"${NSARG[@]}"} delete pod "$POD" --wait=false 2>&1
    ;;

  *) echo "unknown command '$CMD' (up|run|py|down|image)"; exit 2 ;;
esac
