---
name: run_sweep
description: >-
  For a given module (a pipeline stage such as segmentation / transcript_assignment /
  cell_type_annotation / expression_correction, a single method, or 'all'), verify every
  swept method's build is deploy-fresh (via check-component), then launch each method's
  parameter-sweep nebius run, then export a <module>_runs.csv table of (method,
  workflow_id, watch_url). Use when the user wants to kick off a stage's parameter sweeps
  on the cloud, "start all the <stage> runs", (re)launch a method's sweep, or record the
  submitted workflow IDs — after tune-method has already produced the param_sweep files.
---

# run_sweep

Launch a module's parameter-sweep runs *safely* — never launch onto a stale container —
and leave behind a table of what was submitted. It composes two things that already exist
instead of reinventing them:

- **build gate** — [`check-component`](../check-component/SKILL.md)'s `check_component.sh`
  (exit 0 = fully deployed on `origin/main` → `build/main` → ghcr `build_main`).
- **launchers** — `scripts/run_benchmark/param_sweep/run_test_<method>_nebius.sh`
  (produced by [`tune-method`](../tune-method/SKILL.md); each `tw launch`es one method's
  star-sweep and prints a `Workflow <id> submitted` line + watch URL).

This is the automation of what was done by hand for the segmentation stage
(`param_sweep/segmentation_runs.csv`).

## How to run

```bash
.claude/skills/run_sweep/run_sweep.sh <module> [--methods m1,m2,...] [--check-only] [--force] [--tag <tag>]
```

- `<module>` — a pipeline **stage** (`segmentation`, `transcript_assignment`,
  `cell_type_annotation`, `expression_correction`; dashes tolerated), a `methods_*`
  namespace, a **single method** name (`stardist`), or `all`. It is also the basename of
  the exported CSV → `<stage>/<module>_runs.csv` (see Output).
- `--methods` — override the auto-discovered set (e.g. relaunch just the ones that failed).
- `--check-only` — run the build gate and stop; launch nothing.
- `--force` — launch even if a build is stale (default: **abort and launch nothing** if any
  build is stale). Use only when you knowingly accept an old container.
- `--tag` — image tag for the build check (default `build_main`; a `build/<branch>` deploy
  uses `build_<branch>`).

Methods are discovered from the launchers present in `param_sweep/` and filtered to the
requested stage's namespace, so newly-added methods are picked up automatically.

## What it does, in order

1. **Resolve** the module to its method set + CSV basename.
2. **[1/3] Build gate** — `check_component.sh <method>` for each. If any is stale: with
   `--check-only` or without `--force`, print which stage failed and **exit without
   launching**; with `--force`, warn and continue.
3. **[2/3] Launch** — `bash run_test_<method>_nebius.sh` per method, parsing the
   `Workflow <id> submitted` line and the `…/watch/<id>` URL from its output.
4. **[3/3] Export** — write the runs table (`method,workflow_id,watch_url`, successful
   submissions only) to the external results tree, print it, and report any that failed to
   launch (non-zero exit if so).

## Where results go (convention)

All sweep-generated files live **outside the git repo**, under
`~/projects/txsim_results/param_sweep_tests/<pipeline-stage>/` (override the base with
`$TXSIM_RESULTS_DIR`), one folder per stage:

- **run tables** (this skill) → `<stage>/<module>_runs.csv`
- **parsed run results** (the [`parse-sweep-results`](../parse-sweep-results/SKILL.md)
  skill, given a workflow ID) → `<stage>/<method>_run_<id>_results.csv`

A run spanning a single stage (a stage keyword, a namespace, or a single method) writes into
that stage's folder; a multi-stage `all`/`--methods` run writes to the tree root.

## Prerequisites & gotchas

- **Params must be pushed first.** Each launcher reads its `param_sweep/<method>_params.yaml`
  from a **raw-GitHub URL on the current branch** at runtime (see tune-method's params-file
  gotcha). The launcher `curl`-checks this and fails fast; if a launch reports a 404, commit
  **and push** the params yaml (and confirm the URL path includes `param_sweep/`).
- **`tw` must be authenticated** to the Seqera workspace; launches go to whatever
  compute-env/labels each launcher hardcodes (GPU methods keep the `gpu` label + GPU env).
- **Build gate semantics** — the source check compares the component's config on
  `origin/main` to your working tree, so an edit that's only on a feature branch (e.g.
  `fixes`) reads as *not deployed*: the `build/main` pipeline code the run actually uses
  won't contain it. That is correct — merge to `main` and let CI redeploy first. See
  check-component's "Interpreting results" for per-stage fixes.
- The script only records **submitted** runs; it does not wait for or judge their outcome.

## Output

Report to the user: which methods passed/failed the build gate, the workflow IDs + watch
URLs submitted, the CSV path written, and any launch failures (with the likely cause, e.g.
un-pushed params or auth). Don't reprint the whole script.
