---
name: fetch-run-results
description: >-
  Given a Seqera/nebius workflow ID for a task_ist_preprocessing run, copy the run's
  published result files (score_uns.yaml, method_configs.yaml, dataset_uns.yaml, state.yaml,
  ... — the *.yaml the benchmark writes to publish_dir) off the cluster scratch PVC to a
  local folder. Use when the user wants to fetch/download/extract/copy a run's result yamls
  (or its whole output) to a certain folder, pull the scores for a finished run, or get the
  published outputs of a sweep run locally by its workflow ID. Complements parse-sweep-results
  (which reads the task list, not the output files).
---

# fetch-run-results

Copy a finished run's **published output files** from the cluster to a local folder, keyed
off the **workflow ID** that `run_sweep` recorded in `<stage>/<module>_runs.csv`.

A benchmark run does not publish to S3 or anywhere the launch host can see — it writes its
results to `publish_dir` on the shared scratch PVC (`tower-scratch`, mounted at `/scratch`,
namespace `tower-nf`). So there are two steps, both automated here:

1. **Ask Seqera** for the run's `publish_dir` (`tw runs view -i <id> --params`).
2. **`kubectl cp`** the files out of a pod that mounts that PVC (the long-lived `debug`
   pod by default; an ephemeral pod is spawned + torn down if `debug` is gone).

This is the packaged version of the manual copy first done for the cellposev4 run
(`.../segmentation/cellpose4/`).

## What lands in publish_dir

A `run_benchmark` (`-entry auto`) run publishes these top-level yamls (plus a `trace.txt`):

- **`score_uns.yaml`** — the metric **values** per `(dataset, method-combo)`. The one you
  want to rank variants by.
- `method_configs.yaml` / `metric_configs.yaml` — the resolved configs that produced them.
- `dataset_uns.yaml`, `task_info.yaml`, `state.yaml` — dataset/run/output metadata.

By default the skill copies all top-level `*.yaml`; `--all` grabs everything (incl. `trace.txt`).

## How to run

```bash
.claude/skills/fetch-run-results/fetch_results.sh <workflow_id> [<dest_dir>] \
    [--all] [--glob '<pat>'] [--workspace <id>] [--pod <name>] [--keep-pod]
```

- `<workflow_id>` — the run id (e.g. `3JG1eHsWrAeSse`).
- `<dest_dir>` — local folder to copy into (created if absent). The normal case is an
  explicit folder, e.g. `~/projects/txsim_results/param_sweep_tests/segmentation/cellpose4`.
  If omitted, defaults to `<base>/runs/<publish-basename>` (see `--base`).
- `--all` — copy the entire `publish_dir` (files + subdirs, structure preserved), not just
  the top-level `*.yaml`.
- `--glob '<pat>'` — copy top-level files matching this glob instead of `*.yaml`
  (e.g. `--glob 'score_uns.yaml'` for only the scores).
- `--workspace` — Seqera workspace id (default `$TOWER_WORKSPACE_ID` or `167877437119966`).
- `--pod` — pod that mounts the PVC to copy from (default `$SCRATCH_POD` or `debug`).
- `--keep-pod` — if an ephemeral pod was spawned, leave it up (default: delete it).

Example (what the skill was built from):

```bash
.claude/skills/fetch-run-results/fetch_results.sh 3JG1eHsWrAeSse \
    ~/projects/txsim_results/param_sweep_tests/segmentation/cellpose4
```

## How it works

1. **Resolve** `publish_dir` from the run's params (fails fast on auth / wrong id).
2. **Pick a pod** that mounts the scratch PVC: use `--pod`/`debug` if it's Running *and* can
   `test -d publish_dir`; otherwise **spawn** an ephemeral pod (`fetch-results-<id>`) that
   mounts `tower-scratch` using the `process_dataset` image (resolved from `_viash.yaml`, so
   it's a guaranteed-pullable ref), use it, and **tear it down** (unless `--keep-pod`).
3. **Copy** the matching files (or the whole dir with `--all`) into `<dest_dir>` via
   `kubectl cp`, then print the local listing.

## Where results go (convention)

Fits the same external results tree as the sibling skills (override with
`$TXSIM_RESULTS_DIR`), e.g. `~/projects/txsim_results/param_sweep_tests/<stage>/<method>/`:

- **run tables** → `<stage>/<module>_runs.csv` ([`run_sweep`](../run_sweep/SKILL.md))
- **parsed success/fail** → `<stage>/<method>_run_<id>_results.csv`
  ([`parse-sweep-results`](../parse-sweep-results/SKILL.md))
- **fetched output yamls** (this skill) → whatever folder you point it at (e.g.
  `<stage>/<method>/`), or `<base>/runs/<publish-basename>` by default.

## Prerequisites & gotchas

- **`tw` authenticated** to the workspace, and **`kubectl`** pointed at the cluster
  (`kubectl config current-context` → `nebius-mk8s-...`). On Nebius the token expires; if a
  call prints "Switch to your browser to complete authentication", re-auth first.
- **The `debug` pod is the fast path** — it's a long-lived pod in `tower-nf` that already
  mounts `/scratch`, so copies are instant. The ephemeral-pod fallback pulls a multi-GB
  image the first time on an uncached node (minutes), so keep `debug` around if you fetch often.
- **`--all` copies `trace.txt` too** (the Nextflow execution trace, useful for timings/OOM
  diagnosis); the default `*.yaml` leaves it out.
- **"Fetched" ≠ "scored well."** `score_uns.yaml` holds the values; a combo present there
  still may have scored poorly. To turn a run into a ran-vs-failed matrix instead, use
  [`parse-sweep-results`](../parse-sweep-results/SKILL.md).

## Output

Report to the user: the resolved `publish_dir`, which files were copied (and from which pod,
noting if one was spawned/torn down), and the local destination path. Point out that
`score_uns.yaml` is where the metric values live if they want to rank variants next.
