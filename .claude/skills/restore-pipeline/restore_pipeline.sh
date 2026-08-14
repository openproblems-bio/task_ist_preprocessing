#!/usr/bin/env bash
# restore_pipeline.sh — resume a FROZEN or CANCELLED Seqera/Nextflow benchmark run
# on the Nebius k8s env WITHOUT re-running the tasks it already completed.
#
# WHY THIS EXISTS (the core trick):
#   Nextflow reads its -resume cache DB from the head pod's LAUNCH dir (/.nextflow/cache/<session>),
#   NOT from the workDir on /scratch. `tw runs relaunch` DOES pass `-resume <session>` and reuse
#   the workDir, but it starts the new head with an EMPTY /.nextflow/cache and never restores the
#   flushed cache -> so a plain relaunch re-runs EVERYTHING from scratch (verified: 0 Cached).
#   The fix is a `--pre-run` script that copies the real cache into /.nextflow/cache BEFORE
#   Nextflow starts. That + pinning the EXACT original commit (so task hashes match) gives a real
#   resume (verified: 2317/2420 Cached on grave_carson).
#
#   Where the cache lives depends on run state:
#     - RUNNING/FROZEN run  -> ONLY on the head pod's ephemeral fs at /.nextflow/cache/<session>.
#                              If the pod is evicted/cancelled ungracefully it is LOST. Back it up
#                              FIRST (this script does), then cancel.
#     - CANCELLED/finalized -> Tower already flushed it to <workDir>/.nextflow/cache/<session> on
#                              /scratch. Still back it up (insurance) before touching anything.
#
# Usage:
#   restore_pipeline.sh <run_id> [opts]
#
#   --backup-only        Only preserve the cache to /scratch/resume_backup/<run_id>/ and stop.
#                        (Do this IMMEDIATELY for a frozen run, before it is cancelled/evicted.)
#   --cancel             If the run is still live (RUNNING/SUBMITTED), cancel it after backing up
#                        the cache and before relaunching. (In-flight tasks are lost either way —
#                        a frozen head can't cache their completion.)
#   --commit <sha>       Override the revision to relaunch on. Default: the run's own "Commit ID"
#                        (tw runs view). MUST be the exact original commit or every task hash
#                        misses (viash re-stamps meta.git_commit). Passed as --revision (NOT
#                        --commit-id, which does not pin).
#   --session <uuid>     Override the Nextflow session id. Default: the run's "Session ID".
#   --config <file>      Extra Nextflow config layered onto the relaunch (e.g. the fail-fast
#                        retry/memory override — see resume_override.config.template). CACHE-SAFE:
#                        resource/execution directives are not part of the task hash.
#   --params-file <file> Replacement params-file (e.g. to drop a method from a <stage>_methods
#                        list). Build it from `tw runs view -i <run> --params` verbatim, changing
#                        only what you intend. CACHE-SAFE for the branches you keep.
#   --name <name>        Run name for the relaunch (default: <orig-name>_resume).
#   --dry-run            Print the relaunch command (and generated pre-run) without firing it.
#   --workspace <id>     Seqera workspace (default $TOWER_WORKSPACE_ID or 167877437119966).
#   --pod <name>         Preferred pod mounting the scratch PVC (default $SCRATCH_POD or 'debug').
#                        If it's absent, the script AUTO-RESOLVES: it uses any Running pod that
#                        mounts the PVC (prefers an nf-workflow head), else spawns a short-lived
#                        pod that mounts the PVC and tears it down on exit. So 'debug' need not exist.
#   --namespace <ns>     Namespace of the scratch pod / PVC (default $SCRATCH_NS or 'tower-nf').
#   --pvc <name>         Scratch PVC name to mount if a pod must be spawned (default $SCRATCH_PVC
#                        or 'tower-scratch').
#
# Env: TOWER_WORKSPACE_ID, SCRATCH_POD, SCRATCH_NS, SCRATCH_PVC, SCRATCH_MOUNT(/scratch),
#      BACKUP_ROOT(/scratch/resume_backup), RESTORE_IMAGE(ephemeral-pod image)
set -uo pipefail

command -v tw      >/dev/null 2>&1 || { echo "ERROR: 'tw' (Seqera CLI) not found on PATH"; exit 2; }
command -v kubectl >/dev/null 2>&1 || { echo "ERROR: 'kubectl' not found on PATH"; exit 2; }

RUN="" ; BACKUP_ONLY=0 ; DO_CANCEL=0 ; COMMIT="" ; SESSION="" ; CONFIG="" ; PARAMS="" ; NAME="" ; DRY=0
WS="${TOWER_WORKSPACE_ID:-167877437119966}"
POD="${SCRATCH_POD:-debug}"
NS="${SCRATCH_NS:-tower-nf}"
PVC="${SCRATCH_PVC:-tower-scratch}"
MOUNT="${SCRATCH_MOUNT:-/scratch}"
BACKUP_ROOT="${BACKUP_ROOT:-/scratch/resume_backup}"

while [ $# -gt 0 ]; do
  case "$1" in
    --backup-only) BACKUP_ONLY=1;;
    --cancel) DO_CANCEL=1;;
    --commit) COMMIT="$2"; shift;;
    --session) SESSION="$2"; shift;;
    --config) CONFIG="$2"; shift;;
    --params-file) PARAMS="$2"; shift;;
    --name) NAME="$2"; shift;;
    --dry-run) DRY=1;;
    --workspace) WS="$2"; shift;;
    --pod) POD="$2"; shift;;
    --namespace) NS="$2"; shift;;
    --pvc) PVC="$2"; shift;;
    -h|--help) sed -n '2,51p' "$0"; exit 0;;
    -*) echo "unknown option: $1"; exit 2;;
    *) [ -z "$RUN" ] && RUN="$1" || { echo "unexpected arg: $1"; exit 2; };;
  esac
  shift
done
[ -n "$RUN" ] || { echo "ERROR: run_id required. See --help."; exit 2; }
export TOWER_WORKSPACE_ID="$WS"

echo ">>> Resolving run $RUN ..."
VIEW="$(tw runs view -i "$RUN" 2>&1)" || { echo "$VIEW"; echo "ERROR: tw runs view failed (auth? bad id?)"; exit 3; }
get(){ echo "$VIEW" | awk -F'|' "/ $1 /{gsub(/^ +| +\$/,\"\",\$2);print \$2; exit}"; }
STATUS="$(get Status)"; ONAME="$(get 'Run name')"; WORKDIR="$(get Workdir)"
[ -n "$COMMIT" ]  || COMMIT="$(get 'Commit ID')"
[ -n "$SESSION" ] || SESSION="$(get 'Session ID')"
[ -n "$NAME" ]    || NAME="${ONAME}_resume"
echo "    name=$ONAME status=$STATUS"
echo "    session=$SESSION commit=$COMMIT workdir=$WORKDIR"
[ -n "$SESSION" ] && [ -n "$COMMIT" ] || { echo "ERROR: could not resolve session/commit"; exit 3; }

DEST="$BACKUP_ROOT/$RUN/cache"
LC="$(echo "$RUN" | tr 'A-Z' 'a-z')"
HEAD="$(kubectl get pods -n "$NS" --no-headers 2>/dev/null | grep "$LC" | grep workflow | awk '{print $1; exit}')"

# Resolve a pod that mounts the scratch PVC (for reading/writing /scratch). Do NOT
# hard-depend on a 'debug' pod that may not exist: prefer $POD if it's up, else ANY
# Running pod that mounts $PVC (usually a live nf-workflow head), else spawn an
# ephemeral pod that mounts $PVC and tear it down on exit.
SPAWNED=0
find_pvc_pod() {
  kubectl get pods -n "$NS" -o json 2>/dev/null | python3 -c '
import json,sys
pvc=sys.argv[1]
best=None
for p in json.load(sys.stdin).get("items",[]):
    if p.get("status",{}).get("phase")!="Running": continue
    for v in p["spec"].get("volumes",[]):
        if ((v.get("persistentVolumeClaim") or {}).get("claimName"))==pvc:
            nm=p["metadata"]["name"]
            # prefer a workflow head pod (stable), else take the first match
            if nm.startswith("nf-workflow-"): print(nm); sys.exit(0)
            best=best or nm
if best: print(best)
' "$PVC" 2>/dev/null
}
resolve_scratch_pod() {
  if [ "$(kubectl get pod "$POD" -n "$NS" -o jsonpath='{.status.phase}' 2>/dev/null)" = Running ]; then
    SPOD="$POD"; echo "    scratch pod: existing '$SPOD'"; return 0
  fi
  SPOD="$(find_pvc_pod)"
  if [ -n "$SPOD" ]; then echo "    scratch pod: '$SPOD' (mounts $PVC)"; return 0; fi
  # spawn ephemeral — image ref from _viash.yaml (guaranteed-pullable, like fetch-run-results)
  SPOD="restore-scratch-$LC"
  ident() { grep -E "$1" _viash.yaml 2>/dev/null | head -1 | sed -E 's/.*:[[:space:]]*//; s/[[:space:]]*$//'; }
  local ORG PROJ REG IMG
  ORG="$(ident '^organization:')"; PROJ="$(ident '^name:')"; REG="$(ident 'docker_registry:')"; REG="${REG:-ghcr.io}"
  IMG="${RESTORE_IMAGE:-$REG/$ORG/$PROJ/data_processors/process_dataset:build_main}"
  echo "    no scratch pod found — spawning ephemeral '$SPOD' (image $IMG) ..."
  cat <<YAML | kubectl apply -f - >/dev/null 2>&1
apiVersion: v1
kind: Pod
metadata: { name: $SPOD, namespace: $NS }
spec:
  restartPolicy: Never
  containers:
  - { name: main, image: $IMG, command: ["sleep","infinity"], volumeMounts: [{ name: scratch, mountPath: $MOUNT }] }
  volumes:
  - { name: scratch, persistentVolumeClaim: { claimName: $PVC } }
YAML
  SPAWNED=1
  kubectl wait --for=condition=Ready "pod/$SPOD" -n "$NS" --timeout=600s >/dev/null 2>&1 \
    || { echo "ERROR: ephemeral scratch pod not Ready"; exit 4; }
}
cleanup_spod() { [ "${SPAWNED:-0}" -eq 1 ] && kubectl delete pod "$SPOD" -n "$NS" --wait=false >/dev/null 2>&1; }
trap cleanup_spod EXIT

echo ">>> Backing up cache for session $SESSION -> $DEST/$SESSION"
if [ -n "$HEAD" ] && kubectl exec -n "$NS" "$HEAD" -- test -d "/.nextflow/cache/$SESSION" 2>/dev/null; then
  echo "    source: LIVE head pod $HEAD (/.nextflow/cache) — this is the only copy for a running run"
  kubectl exec -n "$NS" "$HEAD" -- sh -c "mkdir -p '$DEST' && cp -a '/.nextflow/cache/$SESSION' '$DEST/' && du -sh '$DEST/$SESSION'" \
    || { echo "ERROR: cache backup from head pod failed"; exit 4; }
else
  resolve_scratch_pod
  if kubectl exec -n "$NS" "$SPOD" -- test -d "$WORKDIR/.nextflow/cache/$SESSION" 2>/dev/null; then
    echo "    source: flushed cache on scratch ($WORKDIR/.nextflow/cache) via pod $SPOD"
    kubectl exec -n "$NS" "$SPOD" -- sh -c "mkdir -p '$DEST' && cp -a '$WORKDIR/.nextflow/cache/$SESSION' '$DEST/' && du -sh '$DEST/$SESSION'" \
      || { echo "ERROR: cache backup from scratch failed"; exit 4; }
  elif kubectl exec -n "$NS" "$SPOD" -- test -d "$DEST/$SESSION" 2>/dev/null; then
    echo "    backup already exists at $DEST/$SESSION (reusing)"
    kubectl exec -n "$NS" "$SPOD" -- du -sh "$DEST/$SESSION" 2>/dev/null
  else
    echo "ERROR: cache DB not found on head pod OR in $WORKDIR/.nextflow/cache OR in $DEST."
    echo "       If the head pod was already evicted/cancelled ungracefully the cache is lost."
    exit 4
  fi
fi

if [ "$BACKUP_ONLY" -eq 1 ]; then
  echo ">>> --backup-only: cache preserved at $DEST/$SESSION. Stopping."
  exit 0
fi

if [ "$DO_CANCEL" -eq 1 ] && ! echo "$STATUS" | grep -qE 'SUCCEEDED|FAILED|CANCELLED|UNKNOWN'; then
  echo ">>> Cancelling live run $RUN (in-flight tasks are unrecoverable for resume anyway) ..."
  [ "$DRY" -eq 1 ] && echo "    (dry-run) tw runs cancel -i $RUN" || tw runs cancel -i "$RUN"
fi

# Generate the pre-run cache-injection script. SOURCED into nf-launcher.sh, so a subshell and
# NO `set -e/-u` (set -u leaks and kills the launcher on unbound NXF_XPACK_LICENSE).
PRERUN="$(mktemp -t prerun_inject_"$RUN"_XXXX.sh)"
cat > "$PRERUN" <<EOF
# auto-generated by restore_pipeline.sh — inject resume cache for session $SESSION
( SESS=$SESSION
  SRC=$DEST/\$SESS
  DST=/.nextflow/cache
  echo "[pre-run] injecting resume cache from \$SRC -> \$DST"
  mkdir -p "\$DST"; rm -rf "\$DST/\$SESS"; cp -a "\$SRC" "\$DST/"
  echo "[pre-run] injected: \$(du -sh "\$DST/\$SESS" | cut -f1)" ) \\
  || echo "[pre-run] injection failed (continuing without resume cache)"
EOF
echo ">>> Generated pre-run injector: $PRERUN"

set -- -i "$RUN" --revision "$COMMIT" --name "$NAME" --pre-run "$PRERUN"
[ -n "$CONFIG" ] && set -- "$@" --config "$CONFIG"
[ -n "$PARAMS" ] && set -- "$@" --params-file "$PARAMS"

echo ">>> Relaunch command:"
echo "    tw runs relaunch $*"
if [ "$DRY" -eq 1 ]; then echo ">>> --dry-run: not firing."; exit 0; fi

OUT="$(tw runs relaunch "$@" 2>&1)"; echo "$OUT"
NEWID="$(echo "$OUT" | grep -oE 'Workflow [A-Za-z0-9]+' | awk '{print $2; exit}')"
echo
echo ">>> Verify cache reuse once the new head starts (want mostly 'Cached process'):"
echo "    H=\$(kubectl get pods --no-headers | grep \$(echo ${NEWID:-<newid>} | tr A-Z a-z) | grep workflow | awk '{print \$1}')"
echo "    kubectl logs \"\$H\" | grep -oE 'Cached process|Submitted process' | sort | uniq -c"
echo "    kubectl logs \"\$H\" | grep -i '\\[pre-run\\]'      # confirm the cache was injected"
