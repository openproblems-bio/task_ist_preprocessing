#!/usr/bin/env bash
# render_qmd.sh — render a Quarto .qmd report against a given dataset-paths YAML on
# a pod on the CURRENT kubectl cluster, using a Viash component's container image so
# the report's spatial/single-cell stack (spatialdata, anndata, scanpy, ...) reads
# the data on the shared scratch PVC natively — no ad-hoc reader installs.
#
# Why a pod: the combined datasets live ONLY on the cluster's scratch PVC
# (tower-scratch, ~/scratch/task_ist_preprocessing/...). A local render can't see
# them. So we mount that PVC into a pod running the process_dataset image (which
# already has spatialdata 0.8 / anndata / scanpy / umap-learn / matplotlib) and add
# just the report-rendering pieces the image lacks: quarto, jupyter/ipykernel, and
# any report-specific python deps (spatialdata_plot by default).
#
# Usage:
#   render_qmd.sh up                              # create pod (scratch mounted) + install render deps
#   render_qmd.sh render <file.qmd> <paths.yaml> [out.html]
#                                                 # stage both into the pod, quarto render, copy HTML back
#   render_qmd.sh down                            # delete the pod
#   render_qmd.sh sh '<bash>'                     # run a bash snippet in the pod (debugging)
#   render_qmd.sh image                           # print the resolved image ref
#
# The report reads its config via $DATASET_CONFIG (falls back to dataset_paths.yaml);
# `render` sets it to the staged YAML's basename and runs quarto with cwd = the qmd's
# dir, so the YAML's relative paths resolve exactly as they do locally. The YAML's
# resources_roots must point at the PVC mount (default /scratch/task_ist_preprocessing).
#
# Env overrides (defaults target this repo's Nebius/Seqera setup):
#   RENDER_POD (render-qmd)   SCRATCH_PVC (tower-scratch)   SCRATCH_NS (tower-nf)
#   SCRATCH_MOUNT (/scratch)  RENDER_COMPONENT (data_processors/process_dataset)
#   RENDER_TAG (build_main)   RENDER_IMAGE (full ref; overrides component/tag)
#   QUARTO_VERSION (1.10.18)  RENDER_PIP (spatialdata_plot)  RENDER_KERNEL (python3)
set -uo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "ERROR: not in a git repo"; exit 2; }
cd "$REPO_ROOT" || exit 2

CMD="${1:-}"
[ -z "$CMD" ] && { grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 2; }

# ---- config (env-overridable) ----
POD="${RENDER_POD:-render-qmd}"
NS="${SCRATCH_NS:-tower-nf}"
PVC="${SCRATCH_PVC:-tower-scratch}"
MOUNT="${SCRATCH_MOUNT:-/scratch}"
QVER="${QUARTO_VERSION:-1.10.18}"
KERNEL="${RENDER_KERNEL:-python3}"
# report-specific python deps the base image lacks (jupyter/ipykernel/quarto are always added)
RENDER_PIP="${RENDER_PIP:-spatialdata_plot}"

# ---- image ref from _viash.yaml (same resolution the cluster pulls with) ----
ident() { grep -E "$1" _viash.yaml | head -1 | sed -E 's/.*:[[:space:]]*//; s/[[:space:]]*$//'; }
ORG="$(ident '^organization:')"; PROJ="$(ident '^name:')"
REG="$(ident 'docker_registry:')"; REG="${REG:-ghcr.io}"
COMPONENT="${RENDER_COMPONENT:-data_processors/process_dataset}"
TAG="${RENDER_TAG:-build_main}"
IMAGE="${RENDER_IMAGE:-$REG/$ORG/$PROJ/$COMPONENT:$TAG}"

CACHE="$MOUNT/.render-cache"          # quarto cached here on the PVC so re-setup skips the 147MB download
WORK="/tmp/render-qmd"                # per-render staging dir inside the pod (kept off the shared PVC)

xp() { kubectl exec -n "$NS" "$POD" -- bash -c "$1"; }   # exec a bash snippet in the pod

case "$CMD" in
  image) echo "$IMAGE" ;;

  up)
    echo "Pod:   $POD   (namespace $NS)"
    echo "Image: $IMAGE"
    echo "Mount: PVC $PVC -> $MOUNT"
    if kubectl get pod "$POD" -n "$NS" >/dev/null 2>&1; then
      echo "pod already exists — (re)installing render deps (idempotent)"
    else
      # kubectl run can't attach a PVC volume -> apply a manifest. `sleep infinity`
      # so the pod never completes (you can't exec a completed pod).
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
      echo "waiting for Ready (image may be multi-GB; cached nodes are instant)..."
      kubectl wait --for=condition=Ready "pod/$POD" -n "$NS" --timeout=600s 2>&1 | tail -1
    fi

    echo "installing render deps (quarto $QVER + jupyter + $KERNEL kernel + pip: $RENDER_PIP)..."
    xp '
      set -e
      QVER="'"$QVER"'"; CACHE="'"$CACHE"'"; KERNEL="'"$KERNEL"'"
      QDIR="$CACHE/quarto-$QVER"
      # quarto: cache the tarball extraction on the PVC, symlink into PATH (idempotent)
      if [ ! -x "$QDIR/bin/quarto" ]; then
        echo "  downloading quarto $QVER ..."
        mkdir -p "$CACHE"
        URL="https://github.com/quarto-dev/quarto-cli/releases/download/v$QVER/quarto-$QVER-linux-amd64.tar.gz"
        curl -sSL "$URL" -o /tmp/quarto.tgz
        tar -xzf /tmp/quarto.tgz -C "$CACHE"
        rm -f /tmp/quarto.tgz
      else
        echo "  quarto $QVER already cached at $QDIR"
      fi
      ln -sf "$QDIR/bin/quarto" /usr/local/bin/quarto
      # quarto jupyter engine + report-specific python deps (pip is idempotent)
      echo "  pip install jupyter ipykernel nbclient '"$RENDER_PIP"' ..."
      pip install --quiet --no-input --root-user-action=ignore jupyter ipykernel nbclient '"$RENDER_PIP"'
      python -m ipykernel install --sys-prefix --name "$KERNEL" --display-name "$KERNEL" >/dev/null 2>&1
      echo "  quarto: $(quarto --version)   kernel: $KERNEL registered"
    '
    ;;

  render)
    QMD="${2:?need a .qmd file as arg 2}"
    YAML="${3:?need a dataset-paths .yaml as arg 3}"
    OUT="${4:-}"
    [ -f "$QMD" ]  || { echo "ERROR: qmd not found: $QMD"; exit 2; }
    [ -f "$YAML" ] || { echo "ERROR: yaml not found: $YAML"; exit 2; }
    kubectl get pod "$POD" -n "$NS" >/dev/null 2>&1 || { echo "ERROR: pod $POD not up — run: $0 up"; exit 2; }

    QBASE="$(basename "$QMD")"; YBASE="$(basename "$YAML")"
    STEM="${QBASE%.*}"
    OUT="${OUT:-./$STEM.html}"

    # Per-dataset work dir keyed by the yaml name, so re-rendering another dataset on
    # the same pod never picks up a previous render's leftovers, and so a resume
    # targets the right render. (Renders on one pod are sequential per caller.)
    RSLUG="$(printf '%s' "$YBASE" | tr -c 'A-Za-z0-9_.-' '_')"
    RWORK="$WORK/$RSLUG"
    RFLAG="$RWORK/render.rc"; RLOG="$RWORK/render.log"; RPID="$RWORK/render.pid"
    RUN="$RWORK/.render_run.sh"; HTML="$RWORK/$STEM.html"

    # Whole-section transcript scatters (5-30M points) run for minutes with no stream
    # output, long enough for the kubectl-exec websocket to reset and kill a synchronous
    # exec mid-render. So run quarto DETACHED (nohup, survives the drop) and poll with
    # short execs. The render is RESUMABLE: if this call is interrupted (e.g. a client
    # timeout), just run the SAME command again — it re-attaches to the running render
    # (or copies the finished HTML) instead of restarting it.
    # done  = finished OK (rc 0 + html present); running = live pid; none = never ran,
    # still running-but-no-pid, OR a PRIOR FAILED run (rc != 0 / no html) -> relaunch fresh.
    _state="$(xp "if [ -f '$RFLAG' ]; then if [ \"\$(cat '$RFLAG' 2>/dev/null)\" = 0 ] && [ -f '$HTML' ]; then echo done; else echo none; fi; elif [ -f '$RPID' ] && kill -0 \$(cat '$RPID') 2>/dev/null; then echo running; else echo none; fi" 2>/dev/null | tr -d '[:space:]')"

    if [ "$_state" = none ]; then
      echo "staging $QBASE + $YBASE into $POD:$RWORK ..."
      xp "rm -rf '$RWORK' && mkdir -p '$RWORK'"
      kubectl cp "$QMD"  "$NS/$POD:$RWORK/$QBASE"  2>&1 | grep -v 'Removing leading' || true
      kubectl cp "$YAML" "$NS/$POD:$RWORK/$YBASE"  2>&1 | grep -v 'Removing leading' || true
      # launcher writes quarto's own rc AFTER it exits, so $RFLAG's presence == done.
      xp "cat > '$RUN' <<'SH'
cd '$RWORK'
DATASET_CONFIG='$YBASE' MPLBACKEND=Agg quarto render '$QBASE' --to html
echo \$? > '$RFLAG'
SH"
      echo "launching detached render (nohup; survives exec/websocket drops) ..."
      xp "cd '$RWORK' && rm -f '$RFLAG' '$HTML' && nohup bash '$RUN' </dev/null >'$RLOG' 2>&1 & echo \$! > '$RPID'"
    elif [ "$_state" = running ]; then
      echo "a render for $YBASE is already in progress in $POD — resuming poll (not restarting) ..."
    else
      echo "render for $YBASE already complete in $POD — copying back ..."
    fi

    if [ "$_state" != done ]; then
      echo "polling (whole-section renders can take several min; interruptible + resumable) ..."
      RC=""; _last=""
      for _i in $(seq 1 300); do          # ceiling; a client timeout just means re-run
        sleep 12
        RC="$(xp "cat '$RFLAG' 2>/dev/null" 2>/dev/null | tr -dc '0-9-')"
        _prog="$(xp "grep -oE 'Cell [0-9]+/[0-9]+|Output created' '$RLOG' 2>/dev/null | tail -n1" 2>/dev/null)"
        [ -n "$_prog" ] && [ "$_prog" != "$_last" ] && { echo "  $_prog"; _last="$_prog"; }
        [ -n "$RC" ] && break
        _alive="$(xp "[ -f '$RPID' ] && kill -0 \$(cat '$RPID') 2>/dev/null && echo y || echo n" 2>/dev/null | tr -d '[:space:]')"
        if [ "$_alive" = n ]; then
          sleep 3
          RC="$(xp "cat '$RFLAG' 2>/dev/null" 2>/dev/null | tr -dc '0-9-')"
          [ -n "$RC" ] && break
          echo "ERROR: render process exited without writing rc (likely OOM/kill). Log tail:"
          xp "tail -n 25 '$RLOG'" 2>/dev/null
          exit 4
        fi
      done
      [ -z "$RC" ] && { echo "still rendering after poll ceiling — resumable: re-run the same command. Log: $POD:$RLOG"; exit 5; }
      if [ "$RC" -ne 0 ]; then
        echo "ERROR: quarto render failed (rc=$RC). Log tail:"
        xp "tail -n 25 '$RLOG'" 2>/dev/null
        exit "$RC"
      fi
    fi

    echo "copying HTML back to $OUT ..."
    kubectl cp "$NS/$POD:$HTML" "$OUT" 2>&1 | grep -v 'Removing leading' || true
    if [ -f "$OUT" ]; then
      echo "OK: $OUT ($(du -h "$OUT" | cut -f1))"
    else
      echo "ERROR: render reported success but $STEM.html not found in pod; check $POD:$RWORK"
      exit 3
    fi
    ;;

  down)
    kubectl delete pod "$POD" -n "$NS" --wait=false 2>&1
    ;;

  sh)
    xp "${2:?need a bash snippet as arg 2}"
    ;;

  *) echo "unknown command '$CMD' (up|render|down|sh|image)"; exit 2 ;;
esac
