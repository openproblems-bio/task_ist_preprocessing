#!/usr/bin/env bash
# run_sweep.sh — for a given module (pipeline stage / method set): verify every
# swept method's build is deploy-fresh, then launch its param-sweep nebius run,
# then export a <module>_runs.csv table of (method, workflow_id, watch_url).
#
# It composes two existing pieces rather than reinventing them:
#   * the build gate  -> .claude/skills/check-component/check_component.sh (exit 0 = OK)
#   * the launchers   -> scripts/run_benchmark/param_sweep/run_test_<method>_nebius.sh
#
# Default behaviour is all-or-nothing: if ANY build is stale, nothing is launched
# (fix the build first, or pass --force). Use --check-only to gate without launching.
#
# Usage:
#   run_sweep.sh <module> [--methods m1,m2,...] [--check-only] [--force] [--tag <tag>]
#
#   <module>   segmentation | transcript_assignment | cell_type_annotation |
#              expression_correction   (a pipeline stage; dashes ok)
#              ... or a methods_* namespace, a single method name, or 'all'.
#              Also the basename of the exported CSV (<module>_runs.csv).
#   --methods  comma/space list overriding the auto-discovered method set.
#   --check-only  run the build checks and stop (launch nothing).
#   --force       launch even if some builds are stale (default: abort if any stale).
#   --tag         image tag passed to check-component (default build_main).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "ERROR: not in a git repo"; exit 2; }
cd "$REPO_ROOT" || exit 2

SWEEP_DIR="scripts/run_benchmark/param_sweep"
CHECK="$REPO_ROOT/.claude/skills/check-component/check_component.sh"
[ -x "$CHECK" ] || CHECK="$SCRIPT_DIR/../check-component/check_component.sh"

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; }

MODULE="" ; METHODS_OVERRIDE="" ; CHECK_ONLY=0 ; FORCE=0 ; TAG="build_main"
while [ $# -gt 0 ]; do
  case "$1" in
    --methods)    METHODS_OVERRIDE="${2:-}"; shift 2;;
    --check-only) CHECK_ONLY=1; shift;;
    --force)      FORCE=1; shift;;
    --tag)        TAG="${2:-}"; shift 2;;
    -h|--help)    usage; exit 0;;
    -*)           echo "unknown option: $1"; usage; exit 2;;
    *)            if [ -z "$MODULE" ]; then MODULE="$1"; else echo "unexpected arg: $1"; exit 2; fi; shift;;
  esac
done
[ -z "$MODULE" ] && { usage; exit 2; }
[ -x "$CHECK" ] || { echo "ERROR: check_component.sh not found/executable at: $CHECK"; exit 2; }

norm="${MODULE//-/_}"

# method name -> its src namespace (e.g. binning -> methods_segmentation)
method_ns() {
  local m="$1" d
  d="$(find src -type d -name "$m" 2>/dev/null | grep -v '/target/' | head -1)"
  [ -z "$d" ] && return 1
  d="${d#src/}"; echo "${d%/$m}"
}

# every method that has a param-sweep nebius launcher
all_methods() {
  local f b
  for f in "$SWEEP_DIR"/run_test_*_nebius.sh; do
    [ -e "$f" ] || continue
    b="$(basename "$f")"; b="${b#run_test_}"; b="${b%_nebius.sh}"
    echo "$b"
  done
}

# ---- resolve method set + CSV basename ----
METHODS=() ; CSV_BASE=""
if [ -n "$METHODS_OVERRIDE" ]; then
  CSV_BASE="${norm#methods_}"
  IFS=', ' read -r -a METHODS <<< "$METHODS_OVERRIDE"
elif [ -f "$SWEEP_DIR/run_test_${MODULE}_nebius.sh" ]; then
  METHODS=("$MODULE"); CSV_BASE="$MODULE"                 # single method
elif [ "$norm" = "all" ]; then
  CSV_BASE="all"; while read -r m; do METHODS+=("$m"); done < <(all_methods)
else
  want_ns="$norm"; [[ "$want_ns" == methods_* ]] || want_ns="methods_$norm"
  CSV_BASE="${want_ns#methods_}"
  while read -r m; do
    ns="$(method_ns "$m" || true)"
    [ "$ns" = "$want_ns" ] && METHODS+=("$m")
  done < <(all_methods)
fi
[ "${#METHODS[@]}" -eq 0 ] && { echo "ERROR: no sweep launchers matched module '$MODULE'"; exit 2; }

# ---- output location: external results tree, organised by pipeline stage ----
# (kept OUT of the git repo; override the base with $TXSIM_RESULTS_DIR)
RESULTS_BASE="${TXSIM_RESULTS_DIR:-$HOME/projects/txsim_results/param_sweep_tests}"
STAGES="$(for m in "${METHODS[@]}"; do ns="$(method_ns "$m" 2>/dev/null || true)"; [ -n "$ns" ] && echo "${ns#methods_}"; done | sort -u)"
n_stages="$(printf '%s\n' "$STAGES" | grep -c .)"
if [ "$n_stages" -eq 1 ]; then RESULTS_DIR="$RESULTS_BASE/$STAGES"; else RESULTS_DIR="$RESULTS_BASE"; fi
mkdir -p "$RESULTS_DIR"
CSV="$RESULTS_DIR/${CSV_BASE}_runs.csv"
echo "Module:  $MODULE   ->  methods: ${METHODS[*]}"
echo "Table:   $CSV"
echo

# ---- [1/3] build gate ----
echo "===== [1/3] build checks (check-component, tag=$TAG) ====="
stale=()
for m in "${METHODS[@]}"; do
  echo "----- $m -----"
  if "$CHECK" "$m" "$TAG"; then
    echo ">> $m: build OK"
  else
    rc=$?
    echo ">> $m: build NOT fresh (check-component rc=$rc)"
    stale+=("$m")
  fi
  echo
done

if [ "${#stale[@]}" -gt 0 ]; then
  echo "Stale / not-deployed: ${stale[*]}"
  [ "$CHECK_ONLY" -eq 1 ] && exit 1
  if [ "$FORCE" -ne 1 ]; then
    echo "ABORTING — launching nothing. Fix the build(s) above, or re-run with --force."
    exit 1
  fi
  echo "WARNING: --force set; launching despite stale build(s)."
  echo
fi
if [ "$CHECK_ONLY" -eq 1 ]; then
  echo "All builds OK (--check-only). Nothing launched."
  exit 0
fi

# ---- [2/3] launch (all builds OK here, or --force) ----
echo "===== [2/3] launching sweeps ====="
rows=() ; failed=() ; launched=0
for m in "${METHODS[@]}"; do
  echo "----- launching $m -----"
  out="$(bash "$SWEEP_DIR/run_test_${m}_nebius.sh" 2>&1)"; rc=$?
  echo "$out"
  wid="$(printf '%s\n' "$out"  | grep -oE 'Workflow[[:space:]]+[A-Za-z0-9]+' | awk '{print $2}' | head -1)"
  wurl="$(printf '%s\n' "$out" | grep -oE 'https://[^[:space:]]+/watch/[A-Za-z0-9]+' | head -1)"
  if [ -n "$wid" ]; then
    rows+=("$m,$wid,$wurl"); launched=$((launched+1))
    echo ">> $m submitted: $wid"
  else
    failed+=("$m")
    echo ">> $m LAUNCH FAILED (rc=$rc, no workflow id in output)"
  fi
  echo
done

# ---- [3/3] export table ----
echo "===== [3/3] exporting table ====="
{
  echo "method,workflow_id,watch_url"
  for r in ${rows[@]+"${rows[@]}"}; do echo "$r"; done
} > "$CSV"
echo "Wrote $CSV"
cat "$CSV"
echo
echo "Launched $launched / ${#METHODS[@]}"
[ "${#failed[@]}" -gt 0 ] && { echo "FAILED to launch: ${failed[*]}"; exit 1; }
exit 0
