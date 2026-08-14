#!/usr/bin/env bash
# sync_scratch_datasets_to_s3.sh — push the combined datasets on the cluster
# scratch PVC up to a NEW `datasets/zarr3/` prefix in S3.
#
#   SRC (cluster):  /scratch/task_ist_preprocessing/datasets/<dataset>/...
#   DST (S3):       s3://openproblems-data/resources/task_ist_preprocessing/datasets/zarr3/<dataset>/...
#
# WHY A POD: the 375 GB of datasets live on the `tower-scratch` PVC (RWX, ns
# tower-nf, mounted at /scratch), which the launch host cannot see; and the AWS
# creds live in your local `op` profile, which the cluster does not have. So we
# bridge them: read the `op` creds locally, hand them to a short-lived pod (via a
# k8s Secret, so they are NOT baked into the pod spec or exec argv), and run
# `aws s3 sync` from inside that pod — data goes cluster -> S3 directly (never
# through your laptop).
#
# ADDITIVE + RESUMABLE: writes to a fresh `zarr3/` prefix and never uses
# --delete, so nothing already in S3 is touched. `aws s3 sync` skips files
# already uploaded with a matching size, so if the transfer drops just re-run
# `start` (or use `status` to watch the one already running).
#
# The sync runs DETACHED inside the pod (nohup -> a log on the PVC), so a dropped
# kubectl connection does not kill a multi-hour transfer.
#
# Usage:
#   scripts/sync_scratch_datasets_to_s3.sh start     [opts]   # spawn pod + launch sync (default)
#   scripts/sync_scratch_datasets_to_s3.sh status    [opts]   # pod state + live sync log + S3 tally
#   scripts/sync_scratch_datasets_to_s3.sh dryrun     [opts]  # aws s3 sync --dryrun (lists, transfers nothing)
#   scripts/sync_scratch_datasets_to_s3.sh cleanup   [opts]   # delete the sync pod + secret
#
#   --only <name>        sync just one dataset dir (repeatable via comma list)
#   --exclude <glob>     pass through to `aws s3 sync --exclude` (repeatable)
#   --wait               after start, block and stream the log until the sync exits
#   --keep               (cleanup) leave the log dir on the PVC (default keeps it anyway)
#   --profile <name>     local AWS profile to read creds from   (default $AWS_PROFILE or 'op')
#   --dst <s3uri>        override destination prefix
#   --concurrency <n>    aws s3 max_concurrent_requests inside the pod (default 64; zarr = many small files)
#   --namespace/--pvc/--mount/--image  cluster knobs (sane defaults below)
set -uo pipefail

ACTION="start"
case "${1:-}" in start|status|dryrun|cleanup) ACTION="$1"; shift;; -*) : ;; "") : ;; *) echo "unknown action: $1"; exit 2;; esac

PROFILE="${AWS_PROFILE:-op}"
NS="${SCRATCH_NS:-tower-nf}"
PVC="${SCRATCH_PVC:-tower-scratch}"
MOUNT="${SCRATCH_MOUNT:-/scratch}"
IMAGE="${SYNC_IMAGE:-amazon/aws-cli:latest}"
SRC_REL="task_ist_preprocessing/datasets"
DST="s3://openproblems-data/resources/task_ist_preprocessing/datasets/zarr3"
POD="${SYNC_POD:-s3sync-datasets}"
SECRET="${POD}-creds"
CONC=64
WAIT=0
ONLY=""
EXCLUDES=()

while [ $# -gt 0 ]; do
  case "$1" in
    --only) ONLY="${ONLY:+$ONLY,}$2"; shift;;
    --exclude) EXCLUDES+=( "$2" ); shift;;
    --wait) WAIT=1;;
    --keep) :;;
    --profile) PROFILE="$2"; shift;;
    --dst) DST="$2"; shift;;
    --concurrency) CONC="$2"; shift;;
    --namespace) NS="$2"; shift;;
    --pvc) PVC="$2"; shift;;
    --mount) MOUNT="$2"; shift;;
    --image) IMAGE="$2"; shift;;
    -h|--help) sed -n '2,50p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 2;;
  esac
  shift
done

SRC="$MOUNT/$SRC_REL"
LOGDIR="$MOUNT/.s3sync"
LOG="$LOGDIR/${POD}.log"
DONE="$LOGDIR/${POD}.exit"

command -v kubectl >/dev/null 2>&1 || { echo "ERROR: kubectl not found"; exit 2; }
command -v aws     >/dev/null 2>&1 || { echo "ERROR: aws CLI not found (needed to read the '$PROFILE' creds)"; exit 2; }
kx(){ kubectl -n "$NS" "$@"; }

# ---------- status ----------
if [ "$ACTION" = "status" ]; then
  echo "=== pod ==="; kx get pod "$POD" -o wide 2>&1 | grep -vE '^$' || true
  if kx get pod "$POD" >/dev/null 2>&1; then
    echo "=== sync log tail ==="; kx exec "$POD" -- sh -c "tail -n 25 '$LOG' 2>/dev/null || echo '(no log yet)'"
    echo "=== exit marker ==="; kx exec "$POD" -- sh -c "[ -f '$DONE' ] && echo \"FINISHED rc=\$(cat '$DONE')\" || echo 'still running (no exit marker)'"
    echo "=== objects under dst so far ==="; kx exec "$POD" -- sh -c "aws s3 ls --recursive '$DST/' 2>/dev/null | wc -l | xargs echo 'S3 objects:'"
  fi
  exit 0
fi

# ---------- cleanup ----------
if [ "$ACTION" = "cleanup" ]; then
  kx delete pod "$POD" --ignore-not-found --wait=false
  kx delete secret "$SECRET" --ignore-not-found
  echo "Deleted pod/$POD and secret/$SECRET (log kept at $LOG on the PVC)."
  exit 0
fi

# ---------- start / dryrun: read creds, ensure secret + pod ----------
AK="$(aws configure get aws_access_key_id     --profile "$PROFILE" 2>/dev/null)"
SK="$(aws configure get aws_secret_access_key --profile "$PROFILE" 2>/dev/null)"
RG="$(aws configure get region                --profile "$PROFILE" 2>/dev/null)"; RG="${RG:-us-west-2}"
[ -n "$AK" ] && [ -n "$SK" ] || { echo "ERROR: could not read aws creds for profile '$PROFILE'"; exit 3; }

echo ">>> Ensuring k8s secret $SECRET (creds from local profile '$PROFILE') ..."
kx delete secret "$SECRET" --ignore-not-found >/dev/null 2>&1
kx create secret generic "$SECRET" \
  --from-literal=AWS_ACCESS_KEY_ID="$AK" \
  --from-literal=AWS_SECRET_ACCESS_KEY="$SK" \
  --from-literal=AWS_DEFAULT_REGION="$RG" >/dev/null || { echo "ERROR: secret create failed"; exit 4; }

if ! kx get pod "$POD" >/dev/null 2>&1; then
  echo ">>> Spawning sync pod $POD (image $IMAGE, mounts $PVC at $MOUNT) ..."
  cat <<YAML | kx apply -f - >/dev/null
apiVersion: v1
kind: Pod
metadata:
  name: $POD
  namespace: $NS
  labels: { app: s3sync-datasets }
spec:
  restartPolicy: Never
  containers:
  - name: sync
    image: $IMAGE
    command: ["/bin/sh","-c","sleep infinity"]
    envFrom:
    - secretRef: { name: $SECRET }
    volumeMounts:
    - name: scratch
      mountPath: $MOUNT
  volumes:
  - name: scratch
    persistentVolumeClaim: { claimName: $PVC }
YAML
  echo ">>> Waiting for pod Ready ..."
  kx wait --for=condition=Ready "pod/$POD" --timeout=300s || { echo "ERROR: pod not Ready"; kx describe pod "$POD" | tail -20; exit 5; }
fi

# preflight: creds + egress + dest reachable (fail fast, cheap)
echo ">>> Preflight: checking S3 write access to $DST ..."
kx exec "$POD" -- sh -c "
  aws configure set default.s3.max_concurrent_requests $CONC
  aws configure set default.s3.max_queue_size 10000
  echo ok | aws s3 cp - '$DST/.synctest' >/dev/null 2>&1 && aws s3 rm '$DST/.synctest' >/dev/null 2>&1 && echo '   write OK' || { echo '   ERROR: cannot write to $DST (creds/egress/bucket policy)'; exit 1; }
" || exit 6

# build the sync command (per-dataset loop for progress; resumable; additive)
EXC=""
for g in "${EXCLUDES[@]:-}"; do [ -n "$g" ] && EXC="$EXC --exclude '$g'"; done
DRY=""; [ "$ACTION" = "dryrun" ] && DRY="--dryrun"

# which datasets
if [ -n "$ONLY" ]; then
  LIST="$(echo "$ONLY" | tr ',' ' ')"
else
  LIST="$(kx exec "$POD" -- sh -c "ls -1 '$SRC' 2>/dev/null")"
fi
[ -n "$LIST" ] || { echo "ERROR: no dataset dirs found under $SRC"; exit 7; }
N=$(echo "$LIST" | wc -w | tr -d ' ')
LIST="$(echo "$LIST" | tr '\n' ' ')"   # space-separate: `for d in $LIST` must be one line in the runner
echo ">>> ${ACTION}: $N dataset dir(s)"
echo "    $SRC/<dataset>  ->  $DST/<dataset>"

# Build the runner as a FILE on the PVC and launch from the file. Do NOT inline it
# into `nohup sh -c "<script>"` — the runner has newlines (the dataset list) and
# $(...)/quotes, and that nested quoting silently broke the detached launch (the
# job never ran yet the script falsely reported success). Feeding the body via
# stdin (`kubectl exec -i ... cat > file`) sidesteps all quoting.
RUNNER="$LOGDIR/${POD}.run.sh"
runner_body() {
  cat <<EOF
mkdir -p '$LOGDIR'
echo "[\$(date -u +%FT%TZ)] START $ACTION $N datasets -> $DST"
rc=0
for d in $LIST; do
  s="$SRC/\$d"; t="$DST/\$d"
  [ -d "\$s" ] || { echo "[skip] \$d (not a dir on scratch)"; continue; }
  echo "[\$(date -u +%FT%TZ)] sync \$d"
  aws s3 sync $DRY "\$s" "\$t" --no-progress $EXC || rc=1
done
echo "[\$(date -u +%FT%TZ)] DONE rc=\$rc"
echo \$rc > '$DONE'
EOF
}
echo ">>> Writing runner to $RUNNER (in-pod) ..."
runner_body | kx exec -i "$POD" -- sh -c "mkdir -p '$LOGDIR'; cat > '$RUNNER'"

if [ "$ACTION" = "dryrun" ]; then
  echo ">>> Running --dryrun in the foreground (no data transferred) ..."
  kx exec "$POD" -- sh "$RUNNER"
  exit 0
fi

# start: launch DETACHED from the file (survives kubectl disconnect via nohup + </dev/null)
echo ">>> Launching detached sync inside the pod (log: $LOG on the PVC) ..."
kx exec "$POD" -- sh -c "rm -f '$DONE' '$LOG'; nohup sh '$RUNNER' > '$LOG' 2>&1 </dev/null & echo \"launched pid \$!\""

# Verify it ACTUALLY started, by LOG GROWTH (the minimal aws-cli image has no
# `ps`, so don't probe processes — an earlier version false-negatived on that).
echo ">>> Verifying the sync is really running ..."
kx exec "$POD" -- sh -c "
  i=0; while [ ! -s '$LOG' ] && [ \$i -lt 10 ]; do i=\$((i+1)); sleep 1; done
  a=\$(wc -c < '$LOG' 2>/dev/null || echo 0); sleep 3; b=\$(wc -c < '$LOG' 2>/dev/null || echo 0)
  if [ -f '$DONE' ]; then echo \"   finished already rc=\$(cat '$DONE')\"; tail -4 '$LOG'
  elif [ \"\$b\" -gt \"\$a\" ]; then echo '   OK — running (log growing). Recent:'; tail -4 '$LOG'
  else echo '   ?? log present but not growing yet (aws may be listing a large dir); check status shortly:'; tail -6 '$LOG'
  fi
"
echo
echo ">>> Sync confirmed running; it survives disconnects. Monitor / finish with:"
echo "    $0 status      # pod state + log tail + S3 object count"
echo "    $0 cleanup     # delete the pod + secret when the log shows DONE rc=0"

if [ "$WAIT" -eq 1 ]; then
  echo ">>> --wait: streaming until exit marker appears ..."
  kx exec "$POD" -- sh -c "until [ -f '$DONE' ]; do sleep 15; done; tail -n 40 '$LOG'; echo; echo \"EXIT rc=\$(cat '$DONE')\""
fi
