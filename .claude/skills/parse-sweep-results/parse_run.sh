#!/usr/bin/env bash
# parse_run.sh — parse a Seqera/nebius sweep run into a per-dataset x combo
# end-to-end success/failure table, saved under the results tree by pipeline stage.
#
# Usage:
#   parse_run.sh <workflow_id> [--workspace <id>] [--stage <stage>] [--method <name>]
#                              [--base <dir>] [--keep-json] [--diagnose [--pod NAME] [--lines N]]
#
#   <workflow_id>  the Seqera run id (e.g. 3JG1eHsWrAeSse), as recorded in <stage>_runs.csv.
#   --workspace    Seqera workspace numeric id (default $TOWER_WORKSPACE_ID or 167877437119966).
#   --stage        override the output stage folder (else inferred from the swept token).
#   --method       override the swept method name used in the filename (else inferred).
#   --base         results-tree base (default $TXSIM_RESULTS_DIR or
#                  ~/projects/txsim_results/param_sweep_tests).
#   --keep-json    keep the raw tasks JSON next to the CSV.
#   --diagnose     pull the tail of the .command.err traceback off the scratch PVC (needs
#                  kubectl) for combos flagged for troubleshooting = methods that FAIL ON THEIR
#                  DEFAULT variant; saves to <method>_run_<id>_failures/*.err.
#   --all-failures with --diagnose, fetch tracebacks for EVERY failed combo, not just the
#                  fail-on-default ones (param-only failures included).
#   --pod NAME     pod (with tower-scratch mounted) to exec for --diagnose (default 'debug').
#   --lines N      lines of .command.err to show/save per failure (default 40).
#
# Output: <base>/<stage>/<method>_run_<id>_results.csv  (+ summary + matrix; on failures also
#         a *_failures.csv with a method-impl-vs-param `likely_class`, and, with --diagnose,
#         the .command.err tracebacks). See SKILL.md "Classifying failures".
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v tw  >/dev/null 2>&1 || { echo "ERROR: 'tw' (Seqera CLI) not found on PATH"; exit 2; }
command -v jq  >/dev/null 2>&1 || { echo "ERROR: 'jq' not found on PATH"; exit 2; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' not found on PATH"; exit 2; }

ID="" ; WS="${TOWER_WORKSPACE_ID:-167877437119966}" ; STAGE="" ; METHOD="" ; KEEP=0
DIAGNOSE=0 ; POD="debug" ; LINES=40 ; ALLFAIL=0
BASE="${TXSIM_RESULTS_DIR:-$HOME/projects/txsim_results/param_sweep_tests}"
while [ $# -gt 0 ]; do
  case "$1" in
    --workspace)    WS="${2:-}"; shift 2;;
    --stage)        STAGE="${2:-}"; shift 2;;
    --method)       METHOD="${2:-}"; shift 2;;
    --base)         BASE="${2:-}"; shift 2;;
    --keep-json)    KEEP=1; shift;;
    --diagnose)     DIAGNOSE=1; shift;;
    --all-failures) ALLFAIL=1; shift;;
    --pod)          POD="${2:-}"; shift 2;;
    --lines)        LINES="${2:-}"; shift 2;;
    -h|--help)   sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    -*)          echo "unknown option: $1"; exit 2;;
    *)           if [ -z "$ID" ]; then ID="$1"; else echo "unexpected arg: $1"; exit 2; fi; shift;;
  esac
done
[ -z "$ID" ] && { echo "usage: parse_run.sh <workflow_id> [opts]"; exit 2; }
export TOWER_WORKSPACE_ID="$WS"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/parse_run.XXXXXX")"
trap '[ "$KEEP" -eq 1 ] || rm -rf "$tmp"' EXIT

# ---- run-level headline ----
status_txt="$(tw runs view -i "$ID" --status 2>&1)" || { echo "$status_txt"; echo "ERROR: could not fetch run $ID (auth? wrong id/workspace?)"; exit 1; }
overall="$(printf '%s\n' "$status_txt" | awk -F'|' '$1 ~ /[[:space:]]Status[[:space:]]/ {gsub(/^[ ]*|[ ]*$/,"",$2); print $2; exit}')"
run_name="$(printf '%s\n' "$status_txt" | awk -F'|' '$1 ~ /Run name/ {gsub(/^[ ]*|[ ]*$/,"",$2); print $2; exit}')"
[ -z "$overall" ] && overall="?"
[ -z "$run_name" ] && run_name="?"
echo "Fetching tasks for run $ID (status: $overall) ..."
case "$overall" in
  SUCCEEDED|FAILED|CANCELLED|UNKNOWN|"?") : ;;
  *) echo "NOTE: run is '$overall' — results may be incomplete (still running).";;
esac

# ---- paginate all tasks (100/page cap) ----
off=0
while :; do
  tw -o json runs view -i "$ID" tasks -c taskId,process,tag,status,exit --max 100 --offset "$off" \
      > "$tmp/page_$off.json" 2>"$tmp/page.err" \
    || { echo "ERROR fetching tasks at offset $off:"; cat "$tmp/page.err"; exit 1; }
  n="$(jq 'length' "$tmp/page_$off.json" 2>/dev/null || echo 0)"
  if [ "$n" -eq 0 ]; then rm -f "$tmp/page_$off.json"; break; fi
  off=$((off + 100))
  [ "$n" -lt 100 ] && break
done
if ls "$tmp"/page_*.json >/dev/null 2>&1; then
  jq -s 'add' "$tmp"/page_*.json > "$tmp/tasks.json"
else
  echo "[]" > "$tmp/tasks.json"
fi
echo "tasks fetched: $(jq 'length' "$tmp/tasks.json")"
echo

# ---- analyse ----
extra=()
[ -n "$STAGE" ]  && extra+=(--stage "$STAGE")
[ -n "$METHOD" ] && extra+=(--method "$METHOD")
py_out="$(python3 "$SCRIPT_DIR/parse_run.py" \
  --tasks "$tmp/tasks.json" --id "$ID" --base "$BASE" \
  --overall "$overall" --run-name "$run_name" \
  ${extra[@]+"${extra[@]}"} 2>&1)"; rc=$?
printf '%s\n' "$py_out"

OUT_DIR="$(printf  '%s\n' "$py_out" | sed -n 's/^OUT_DIR=//p'      | tail -1)"
FT_FILE="$(printf  '%s\n' "$py_out" | sed -n 's/^FAILED_TASKS=//p' | tail -1)"
CSV_FILE="$(printf '%s\n' "$py_out" | sed -n 's/^WROTE_CSV=//p'    | tail -1)"

# ---- --diagnose: pull .command.err tracebacks off scratch for flagged (fail-on-default)
#      combos (or all failures with --all-failures) ----
if [ "$DIAGNOSE" -eq 1 ] && [ -n "$FT_FILE" ] && [ -f "$FT_FILE" ]; then
  if ! command -v kubectl >/dev/null 2>&1; then
    echo; echo "WARN: --diagnose needs kubectl (not found) — skipping traceback fetch."
  else
    base="$(basename "${CSV_FILE:-run_${ID}_results.csv}" _results.csv)"
    errdir="${OUT_DIR:-.}/${base}_failures"
    mkdir -p "$errdir"
    scope="fail-on-default only"; [ "$ALLFAIL" -eq 1 ] && scope="all failures"
    echo; echo "=== DIAGNOSE ($scope): .command.err tail via pod '$POD' ==="
    while IFS=$'\t' read -r tid ds combo step ex flag; do
      [ -z "$tid" ] && continue
      [ "$ALLFAIL" -eq 1 ] || [ "${flag:-0}" = "1" ] || continue   # default: flagged only
      wd="$(tw runs view -i "$ID" task -t "$tid" 2>/dev/null | awk -F'|' '/Work directory/{gsub(/^[ ]*|[ ]*$/,"",$2);print $2}')"
      safe="$(printf '%s' "${ds}__${combo}" | tr ' /|=' '_____')"
      errf="$errdir/${safe}.err"
      echo; echo "----- $ds | $combo   (step=$step exit=$ex taskId=$tid) -----"
      if [ -z "$wd" ]; then echo "  (no workdir from tw)"; continue; fi
      if kubectl exec "$POD" -- sh -c "tail -n $LINES '$wd/.command.err' 2>/dev/null" > "$errf" 2>/dev/null && [ -s "$errf" ]; then
        sed 's/^/  /' "$errf"; echo "  (saved: $errf)"
      else
        echo "  (could not read $wd/.command.err via pod '$POD')"
      fi
    done < "$FT_FILE"
  fi
fi

[ -n "$FT_FILE" ] && rm -f "$FT_FILE" 2>/dev/null   # transient hand-off file
[ "$KEEP" -eq 1 ] && echo "(--keep-json) raw tasks JSON kept at: $tmp/tasks.json"
exit $rc
