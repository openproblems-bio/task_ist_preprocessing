---
name: parse-sweep-results
description: >-
  Given a Seqera/nebius workflow ID for a task_ist_preprocessing parameter-sweep run,
  parse its task list into a per-dataset x method-combo end-to-end success/failure table
  (which param combos and which methods succeeded vs failed through to the final metric),
  classify each failure (method-implementation error vs param-combo error) and FLAG for
  troubleshooting only the methods that fail on their DEFAULT variant, optionally pulling the
  real .command.err tracebacks. Saves under the results tree as <stage>/<method>_run_<id>_*.csv.
  Use when the user wants to parse/analyse a sweep run's results by its workflow ID, see which
  combinations succeeded or failed, or triage/troubleshoot the failures. Complements run_sweep
  (which launches the runs and records the IDs).
---

# parse-sweep-results

Turn a finished (or running) sweep into a results table, keyed off the **workflow ID** that
`run_sweep` recorded in `<stage>/<module>_runs.csv`. It pulls the run's tasks from Seqera
with `tw`, decides for every `(dataset, method-combo)` whether it **succeeded end-to-end**
(reached a **COMPLETED terminal metric** task `…/metric_<name>`), and triages the failures.

## How to run

```bash
.claude/skills/parse-sweep-results/parse_run.sh <workflow_id> \
    [--workspace <id>] [--stage <stage>] [--method <name>] [--base <dir>] [--keep-json] \
    [--diagnose [--all-failures] [--pod NAME] [--lines N]]
```

- `<workflow_id>` — the run id (e.g. `3JG1eHsWrAeSse`).
- `--workspace` — Seqera workspace numeric id (default `$TOWER_WORKSPACE_ID` or `167877437119966`).
- `--stage` / `--method` — override the auto-inferred output folder / filename (normally
  inferred from the swept token, e.g. `segm_cellposev4_*` → stage `segmentation`, method `cellposev4`).
- `--base` — results-tree base (default `$TXSIM_RESULTS_DIR` or `~/projects/txsim_results/param_sweep_tests`).
- `--diagnose` — for the **flagged** failures (methods that fail on their default variant), pull
  the tail of each `.command.err` traceback off the scratch PVC via `kubectl exec` (default pod
  `debug`, which has `tower-scratch` mounted). Saves to `<method>_run_<id>_failures/*.err`.
- `--all-failures` — with `--diagnose`, fetch tracebacks for **every** failed combo (incl. param-only).
- `--pod NAME` / `--lines N` — pod to exec (default `debug`) / lines of traceback (default 40).

## How it works

Task **tags** encode the whole pipeline config `"<dataset>/<seg_token>/<ass_token>/…/metric_<name>"`
(control = `"<dataset>/control_<name>/metric_<name>"`). The script:

1. Fetches the run status headline, then **paginates all tasks** (`tw` caps at 100/page).
2. **Auto-detects the swept axis** = the chain-token position whose value varies most across
   non-control tasks — works for any stage, for a single method's param sweep OR several methods
   at one stage (each method → `param=default`).
3. **Parses each varying token** into `method / param / value` (base-method = the token-prefix
   shared by the most variants; a different method used as the fixed baseline is its own method).
4. Marks each `(dataset, combo)` **SUCCEEDED** iff it has a COMPLETED `metric_*` task, else **FAILED**.
5. **Triages + classifies** every failure (see below), writes a `*_failures.csv`, and with
   `--diagnose` fetches the tracebacks.

### Truncation gotcha (important)

Seqera **truncates long task tags**, so the trailing `/metric_<name>` is chopped for long combos
(worst on long dataset names). Detect terminal-metric tasks by **process name** (discovered from
the short tags that still show a `metric_` segment), never by tag suffix. The swept token near the
front of the tag survives, so combo identity and success are still recovered correctly.

## Classifying failures (triage)

The point of the failure analysis is to separate **fix-the-code** from **bad-input-param**:

- **🚩 TROUBLESHOOT = a method that fails on its DEFAULT variant.** This is the primary triage
  rule: if the method breaks out-of-the-box on a dataset, that's a real implementation/robustness
  problem worth chasing. Methods whose default passes and only a *non-default param* fails are
  **not flagged** (the sweep is *supposed* to probe bad params) — they're recorded as param-only.
  In a multi-method run (no params) every method's sole variant IS its default, so any failure flags it.

Each failed combo also gets a heuristic `likely_class` to explain the mechanism:

- **where it dies:** FAILED in the combo's **own method step** → suspect the method/impl; FAILED in
  a **downstream** step (the method produced output, a later stage choked) → a param yielded
  degenerate output (e.g. empty segmentation → downstream divide-by-zero).
- **failure ratio:** if ~all of a method's variants fail on a dataset → *method/dataset-level* (impl
  bug or data-incompat, not one param); if only this variant fails → *param-specific*.

**Deploy-state caveat:** an "impl error" verdict may just be an **undeployed fix** — the run uses the
`build/main` container, which can lag your working tree. Before blaming the code, cross-check with the
`check-component` skill (this is real: the stardist `min_overlap` fix existed locally but the run ran
the old container). `run_sweep`'s build gate exists precisely to avoid launching onto a stale build.

**Retry semantics** (`src/base/labels_nebius.config`): exit-1 (code error) is retried once (2 attempts),
OOM (137/139) up to 3, then `ignore` so the run completes. A **FAILED** combo therefore failed on
**every** attempt — not transient.

## Output (under `<base>/<stage>/`)

- `<method>_run_<id>_results.csv` — `dataset,combo,stage,method,param,value,reached_metric,end_to_end_status`.
- `<method>_run_<id>_failures.csv` (only if there are failures) —
  `dataset,combo,method,param,value,is_default,failing_step,exit,task_id,troubleshoot,likely_class`
  (`troubleshoot=yes` ⇔ the method fails on its default variant).
- `<method>_run_<id>_failures/*.err` (with `--diagnose`) — the `.command.err` tail per fetched failure.

The console prints the run headline, a combos×datasets PASS/FAIL matrix, and a FAILURE DIAGNOSIS
listing the flagged (`[FLAG]`) methods and each failure's step/exit/`likely_class`.

Report to the user: overall status, per-dataset pass/fail counts, **which methods are flagged for
troubleshooting (fail on default)** vs the param-only failures, and the CSV paths.

## Notes & caveats

- **"Succeeded" = ran to completion, not "scored well."** To rank passing combos by metric *value*,
  use the `fetch-run-results` skill to pull `score_uns.yaml` (offer as a follow-up).
- Multi-axis (hand-built) runs fall back to labelling combos by the full varying-token tuple.
- Requires `tw` (authenticated), `jq`, `python3`; `--diagnose` also needs `kubectl` + a pod with
  `tower-scratch` mounted.
