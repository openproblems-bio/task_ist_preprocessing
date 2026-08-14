---
name: restore-pipeline
description: >-
  Resume a FROZEN or CANCELLED Seqera/Nextflow benchmark run on the Nebius k8s env WITHOUT
  re-running the tasks it already completed — by preserving its -resume cache and injecting it
  into a fresh head via a --pre-run script (Tower's own relaunch does NOT restore the cache, so
  a plain relaunch re-runs everything from scratch). Use when a run's head froze (JVM hung —
  still "RUNNING" but silent for hours, submitting nothing) or was cancelled, and you want to
  recover it; also for relaunching a run with a method dropped or a changed retry/memory policy
  while keeping the cache. Complements fetch-run-results (which copies a finished run's outputs).
---

# restore-pipeline

Resume a run that **froze** or was **cancelled**, reusing everything it already computed. On this
Nebius compute env a plain `tw runs relaunch --resume` re-runs **from scratch** — this skill makes
it a real resume (validated: **2317/2420 tasks Cached** on `grave_carson`, `~1031/1061` on
`nauseous_celsius`).

Background + full war story: memory **[[nextflow-resume-cache-on-ephemeral-head-pod]]** and
**[[tw-resume-git-commit-cache-bust]]**.

## When to use

- A run's **head froze**: Tower still shows `RUNNING`, but the head pod's log is silent for hours,
  `RESTARTS 0`, and it has **submitted no new task pods** since the freeze (a hung JVM — often GC
  death-spiral under `0/76 nodes: Insufficient memory/ephemeral-storage`). It will never finish.
- A run was **cancelled** and you want to continue it.
- You want to **relaunch with a change** (drop a buggy method, change retry/memory) but keep the
  cache.

## The core problem (why a plain relaunch fails)

Nextflow reads its `-resume` cache DB (`.nextflow/cache/<session>/` = LevelDB `db/` +
`index.<runName>`) from the **head pod's LAUNCH dir** (`/` -> `/.nextflow/cache`), **not** from the
workDir on `/scratch`. `tw runs relaunch` passes `-resume <session>` and reuses the workDir, but
starts the new head with an **empty** `/.nextflow/cache` and never restores the flushed cache →
`0 Cached / N Submitted` (verified). The fix: a `--pre-run` script copies the real cache into
`/.nextflow/cache` **before** Nextflow starts.

**Where the cache lives depends on state** (the script handles both):
- **RUNNING / FROZEN** → only on the head pod's **ephemeral** fs (`/.nextflow/cache/<session>`).
  If the pod is evicted or cancelled ungracefully, it is **LOST**. Back it up FIRST.
- **CANCELLED / finalized** → Tower flushed it to `<workDir>/.nextflow/cache/<session>` on
  `/scratch`. Still back it up (insurance).

## How to run

```bash
.claude/skills/restore-pipeline/restore_pipeline.sh <run_id> [opts]
```

**Frozen run — do the urgent backup first, then resume:**
```bash
# 1) IMMEDIATELY preserve the cache off the live (frozen) head pod, before anything kills it:
.claude/skills/restore-pipeline/restore_pipeline.sh <run_id> --backup-only
# 2) resume (cancels the frozen run, injects the cache, pins the exact commit):
.claude/skills/restore-pipeline/restore_pipeline.sh <run_id> --cancel
```

**Cancelled run — one shot:**
```bash
.claude/skills/restore-pipeline/restore_pipeline.sh <run_id>
```

**With modifications (both cache-safe — see below):**
```bash
.claude/skills/restore-pipeline/restore_pipeline.sh <run_id> --cancel \
    --config      resume_override.config \      # fail-fast retry/memory (template provided)
    --params-file my_params.yaml                # e.g. a method dropped
```

Key options: `--backup-only`, `--cancel`, `--commit <sha>`, `--session <uuid>`, `--config <file>`,
`--params-file <file>`, `--name`, `--dry-run`, `--workspace`, `--pod`. Defaults for `--commit` /
`--session` come from `tw runs view -i <run>` (the run's own Commit ID / Session ID).

## What the script does

1. Resolve `session`, `commit`, `status`, `workDir` from `tw runs view`.
2. **Back up the cache** to `/scratch/resume_backup/<run_id>/cache/<session>` — from the live head
   pod if running, else from the flushed copy on `/scratch`.
3. (`--cancel`) cancel the live run.
4. Generate a `--pre-run` injector (subshell; **no `set -e/-u`** — `set -u` leaks into
   `nf-launcher.sh` and kills it on unbound `NXF_XPACK_LICENSE`).
5. `tw runs relaunch -i <run> --revision <commit> --pre-run <injector> [--config …] [--params-file …]`.
6. Print the verification commands.

## Two cache-safe modifications you can fold into the resume

Nextflow hashes a task's *identity* (script, inputs, container, params it consumes, `git_commit`),
**not** resource/execution directives or unrelated branches — so these don't bust the cache:

- **Drop a method** → build a replacement `--params-file` from `tw runs view -i <run> --params`
  *verbatim*, removing the method from its `<stage>_methods` array (the method lists live inside a
  JSON-string `settings` key; the stage keys are real workflow params — `--transcript_assignment_methods`
  etc.). The removed branch generates no tasks; kept branches stay Cached. (Include-list only:
  `checkItemAllowed(item, include=state[key], exclude=null)` — to drop one, list all-except-it.)
- **Change retry/memory** → `--config resume_override.config` (copy the provided
  `resume_override.config.template`). Fail-fast (`errorStrategy='ignore'; maxRetries=0`) stops OOM
  tasks escalating to 480 GB and hogging nodes. Mind the two footguns documented in the template
  (`get_memory` forcing 480 GB when `attempt==maxRetries`; `exitStrat` halting on too-few retries).

## Verify (always)

Once the new head schedules and starts:
```bash
H=$(kubectl get pods --no-headers | grep <newid-lowercase> | grep workflow | awk '{print $1}')
kubectl logs "$H" | grep -i '\[pre-run\]'                                  # cache was injected
kubectl logs "$H" | grep -oE 'Cached process|Submitted process' | sort | uniq -c   # want mostly Cached
```
Mostly `Cached process` → success. `0 Cached` → the injection or the commit pin failed (see gotchas).

## Gotchas

- **Pin the EXACT original commit.** `--revision <sha>` (the script defaults to the run's own
  Commit ID). A newer `build/main` re-stamps `meta.git_commit` into every module → 100% hash miss.
  `--commit-id` does **not** pin — pass the SHA as the revision. See [[tw-resume-git-commit-cache-bust]].
- **In-flight tasks are unrecoverable.** A frozen head can't record a task's completion in the
  cache, so any task still running at freeze time re-runs on resume (its on-disk output isn't in
  the cache DB). Waiting for them is pointless — resume now under a live head that *will* record them.
- **The pre-run must not `set -e/-u`.** It is sourced into `nf-launcher.sh`; `set -u` kills it on
  unbound `NXF_XPACK_LICENSE`. The generated script uses a subshell and never fails the launcher.
- **Capacity.** These runs have huge per-task memory (up to 480 GiB ≈ a whole node). Don't resume
  many at once, and consider the fail-fast `--config` to cap memory (300 GiB) so tails schedule.
- **Auth/context.** `tw` authenticated to the workspace and `kubectl` on the nebius context
  (`kubectl config current-context` → `nebius-mk8s-…`). Nebius tokens expire; re-auth if a call
  prints "Switch to your browser to complete authentication".
- **macOS bash 3.2** has no `declare -A` — any helper loops you write around this must use temp
  files or plain vars, not associative arrays.

## Output

Report: the resolved session/commit, where the cache was backed up (and whether from the live head
or scratch — critical for a frozen run), the new run id, and the Cached/Submitted verification.
Then hold or proceed per the user's plan (e.g. stagger multiple resumes for capacity).
