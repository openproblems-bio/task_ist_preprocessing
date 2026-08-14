#!/usr/bin/env bash
# check_component.sh — verify a Viash component's change is fully deployed.
#
# Runs three checks and prints a verdict:
#   [1/3] source: is the component's config committed to origin/main?
#   [2/3] codegen: did build/main redeploy the latest origin/main commit?
#   [3/3] container: was the ghcr :build_main image rebuilt from the latest commit?
#
# The container check reads the image's org.opencontainers.image.revision label
# straight from ghcr (anonymous pull token) and compares it to origin/main HEAD,
# so it does NOT pull the (multi-GB) image.
#
# Usage:
#   check_component.sh <component>                         # e.g. segger
#   check_component.sh <namespace>/<component>             # e.g. methods_transcript_assignment/segger
#   check_component.sh <path/to/config.vsh.yaml>
#   check_component.sh <component> <image-tag>             # tag defaults to build_main
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "ERROR: not in a git repo"; exit 2; }
cd "$REPO_ROOT" || exit 2

ARG="${1:-}"
TAG="${2:-build_main}"
[ -z "$ARG" ] && { echo "usage: check_component.sh <component|namespace/component|config-path> [image-tag]"; exit 2; }

# ---- resolve the component config path ----
if [[ "$ARG" == *config.vsh.yaml && -f "$ARG" ]]; then
  CFG="$ARG"
elif [[ -f "src/$ARG/config.vsh.yaml" ]]; then
  CFG="src/$ARG/config.vsh.yaml"
elif [[ -d "$ARG" && -f "$ARG/config.vsh.yaml" ]]; then
  CFG="$ARG/config.vsh.yaml"
else
  HITS="$(find src -type d -name "$ARG" 2>/dev/null | while read -r d; do [ -f "$d/config.vsh.yaml" ] && echo "$d"; done)"
  N="$(printf '%s\n' "$HITS" | grep -c . || true)"
  if [ "$N" -eq 0 ]; then echo "ERROR: no component config found for '$ARG'"; exit 2; fi
  if [ "$N" -gt 1 ]; then echo "ERROR: '$ARG' is ambiguous, matches:"; printf '%s\n' "$HITS"; exit 2; fi
  CFG="$HITS/config.vsh.yaml"
fi

COMP_DIR="$(dirname "$CFG")"          # src/<namespace>/<component>
COMP="$(basename "$COMP_DIR")"        # <component>
NS="${COMP_DIR#src/}"; NS="${NS%/$COMP}"   # <namespace> (may be nested)

# ---- project identity from _viash.yaml ----
ident() { grep -E "$1" _viash.yaml | head -1 | sed -E 's/.*:[[:space:]]*//; s/[[:space:]]*$//'; }
ORG="$(ident '^organization:')"
PROJ="$(ident '^name:')"
REG="$(ident 'docker_registry:')"; REG="${REG:-ghcr.io}"
IMG="$ORG/$PROJ/$NS/$COMP"            # repo path after the registry

echo "Component: $COMP  ($COMP_DIR)"
echo "Image:     $REG/$IMG:$TAG"
echo

git fetch origin --quiet 2>/dev/null || echo "(warning: git fetch failed — results may be stale)"
MAIN_SHA="$(git rev-parse origin/main 2>/dev/null || echo '')"
MAIN_SUBJ="$(git log -1 --format=%s origin/main 2>/dev/null || echo '')"

ok_src=1 ok_deploy=1 ok_cont=1

# ================= [1/3] source committed to main =================
echo "[1/3] Source committed to origin/main"
echo "  origin/main HEAD:            ${MAIN_SHA:0:8}  $MAIN_SUBJ"
UNCOMMITTED="$(git status --porcelain -- "$CFG" 2>/dev/null)"
if [ -n "$UNCOMMITTED" ]; then
  echo "  local uncommitted change:    YES  (working-tree edits to config not committed)"
  ok_src=0
else
  echo "  local uncommitted change:    no"
fi
if git diff --quiet origin/main -- "$CFG" 2>/dev/null; then
  echo "  config matches origin/main:  yes"
else
  echo "  config matches origin/main:  NO  (your version differs from what's on main)"
  ok_src=0
fi
[ "$ok_src" -eq 1 ] && echo "  => OK: this component's config on origin/main matches your working tree" \
                     || echo "  => NOT on main: merge/commit the change before it can build"
echo

# ================= [2/3] build/main regenerated =================
echo "[2/3] build/main regenerated (codegen)"
DEPLOY_MSG="$(git log -1 --format=%s origin/build/main 2>/dev/null || echo '')"
DEPLOY_SHA="$(printf '%s' "$DEPLOY_MSG" | grep -oE '[0-9a-f]{40}' | head -1)"
echo "  build/main last deploy:      ${DEPLOY_SHA:0:8}  ($DEPLOY_MSG)"
if [ -n "$DEPLOY_SHA" ] && [ "$DEPLOY_SHA" = "$MAIN_SHA" ]; then
  echo "  matches origin/main HEAD:    yes"
  echo "  => OK: build/main codegen is current"
else
  echo "  matches origin/main HEAD:    NO"
  echo "  => build/main has not redeployed the latest main yet (CI Build may be running)"
  ok_deploy=0
fi
echo

# ================= [3/3] container rebuilt on ghcr =================
echo "[3/3] Container rebuilt on $REG ($TAG)"
CONT_OUT="$(python3 "$SCRIPT_DIR/ghcr_inspect.py" "$REG" "$IMG" "$TAG" "$MAIN_SHA")"
C_STATUS="$(printf '%s\n' "$CONT_OUT" | awk '/^STATUS/{print $2}')"
if [ "$C_STATUS" != "ok" ]; then
  C_DETAIL="$(printf '%s\n' "$CONT_OUT" | awk '/^DETAIL/{print $2}')"
  echo "  image on registry:           NOT FOUND ($C_DETAIL)"
  echo "  => no $TAG image has been pushed for this component"
  ok_cont=0
else
  C_REV="$(printf '%s\n' "$CONT_OUT"  | awk '/^REVISION/{print $2}')"
  C_CREATED="$(printf '%s\n' "$CONT_OUT" | awk '/^CREATED /{print $2}')"
  C_CREATED_EPOCH="$(printf '%s\n' "$CONT_OUT" | awk '/^CREATED_EPOCH/{print $2}')"
  C_BASE="$(printf '%s\n' "$CONT_OUT"  | awk '/^BASE_CREATED/{print $2}')"
  C_LAYERS="$(printf '%s\n' "$CONT_OUT" | awk '/^LAYERS/{print $2}')"
  C_MATCH="$(printf '%s\n' "$CONT_OUT" | awk '/^MATCH/{print $2}')"
  echo "  image revision label:        ${C_REV:0:8}"
  echo "  matches origin/main HEAD:    $([ "$C_MATCH" = yes ] && echo yes || echo 'NO  <-- STALE')"
  echo "  image built:                 $C_CREATED"
  echo "  base layer created:          $C_BASE"
  echo "  layer count:                 $C_LAYERS"
  if [ "$C_MATCH" = yes ]; then
    echo "  => OK: container was built from the current origin/main commit"
  else
    echo "  => STALE: container was built from ${C_REV:0:8}, not ${MAIN_SHA:0:8}."
    echo "     The container push for the latest commit has not succeeded"
    echo "     (CI Build still running, or the image build failed)."
    ok_cont=0
  fi

  # Cross-check: an honestly-built image is CREATED at (or after) the commit it
  # claims to be built from. If 'created' predates the revision's commit time, the
  # layers were NOT rebuilt for that revision — the label was re-stamped over old
  # layers (Docker cache reuse, a re-pushed tag, or a heavy build that failed and
  # left the previous image). The revision label alone can't catch this; the
  # timestamp does. GRACE absorbs build-vs-git clock skew.
  GRACE=120
  REV_EPOCH="$(git show -s --format=%ct "$C_REV" 2>/dev/null)"
  # Show BOTH times in UTC. The image 'created' is UTC (…Z) but git's default
  # commit time carries the author's zone — comparing them as-is is a classic
  # mistake (an image looks "older" than its commit only because of the offset).
  REV_UTC="$(TZ=UTC git show -s --date=format-local:'%Y-%m-%dT%H:%M:%SZ' --format=%cd "$C_REV" 2>/dev/null)"
  if [ -n "$REV_EPOCH" ] && printf '%s' "$C_CREATED_EPOCH" | grep -qE '^[0-9]+$'; then
    echo "  revision commit time (UTC):  $REV_UTC"
    if [ "$C_CREATED_EPOCH" -lt "$((REV_EPOCH - GRACE))" ]; then
      echo "  created-vs-commit:           STALE  <-- image built BEFORE its own revision commit"
      echo "  => STALE: image 'created' ($C_CREATED, UTC) predates the commit ${C_REV:0:8}"
      echo "     it is labelled with ($REV_UTC). The revision label was stamped but the"
      echo "     layers were NOT rebuilt (cache reuse / re-stamped tag / failed heavy build)."
      echo "     Rebuild+push with a clean cache."
      ok_cont=0
    else
      echo "  created-vs-commit:           ok (image built after its revision commit)"
    fi
  else
    echo "  created-vs-commit:           (skipped — revision commit not in local history)"
  fi
fi
echo

# ================= verdict =================
if [ "$ok_src" -eq 1 ] && [ "$ok_deploy" -eq 1 ] && [ "$ok_cont" -eq 1 ]; then
  echo "VERDICT: OK — fully deployed. A run with --pull-latest will use the current build."
  exit 0
else
  echo "VERDICT: NOT fully deployed — a run now would use a stale build. See the failing stage above."
  [ "$ok_src"    -eq 0 ] && echo "  - source not on main: commit/merge the change to main"
  [ "$ok_deploy" -eq 0 ] && echo "  - build/main behind: wait for / re-trigger the CI Build workflow"
  [ "$ok_cont"   -eq 0 ] && echo "  - container stale: check the Build run for this component; if it failed, rebuild+push manually"
  exit 1
fi
