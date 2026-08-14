---
name: tune-method
description: >-
  Study a benchmark method component in depth, identify the parameters that matter for
  optimization (with explicit quality-vs-speed criteria), and set up a parameter-sweep test
  run for it. Three phases: (1) collect knowledge — git history, config-linked docs, upstream
  repo/tutorial, and the ACTUAL publication; (2) identify optimization parameters + ranges
  with criteria, exposing high-value knobs the component doesn't yet surface; (3) generate a
  cloud (nebius) run script + committed params yaml under scripts/run_benchmark/param_sweep/.
  Use when the user wants to
  optimize/tune a method, run a parameter sweep, "study everything on <method>", find which
  knobs matter, or prep a benchmark run that varies a method's parameters.
---

# tune-method

Turning a benchmark method into a *tunable* one is three phases that build on each other.
Do them in order — the parameter choices in phase 2 must be grounded in phase 1's evidence,
and the sweep in phase 3 is only as good as phase 2's criteria.

| Phase | Goal | Output |
|-------|------|--------|
| 1. Collect knowledge | understand the method, its defaults, and the real upstream behaviour | an accurate `NOTES.md` (+ citation caveats) |
| 2. Identify parameters + criteria | decide which knobs matter and over what range, expose missing ones | an "Optimization / tuning" section; possibly new config args |
| 3. Set up the sweep | make it runnable | a cloud `run_test_<method>_nebius.sh` + committed `<method>_params.yaml`, both under `scripts/run_benchmark/param_sweep/` |

Scale the effort to the ask: a quick "what are the knobs" stops after phase 2; "set up a
sweep to optimize X" runs all three.

---

## Phase 1 — Collect knowledge

Resolve the method to `src/<namespace>/<component>/` (e.g. `methods_segmentation/cellposev4`).
Then read, in this order, and **don't skip the paper**:

1. **`config.vsh.yaml`** — the exposed `arguments:` and their defaults, the `engines[docker]`
   setup (version pins — often load-bearing), `links.documentation` / `links.repository`,
   and `references.doi`.
2. **`script.py` / `script.R`** — which args are *actually* forwarded to the underlying
   tool vs hardcoded/ignored, and what the input contract is (e.g. only `image[0]` — a single
   channel — is segmented). The gap between "config arg" and "what the tool call receives" is
   where real behaviour lives.
3. **Git history** — `git log --oneline -- <dir>` and read the introducing commit(s). This is
   how you learn *why* defaults are the way they are — frequently they are **deliberately
   speed-tuned** to fit a CI/GPU time budget (that context is the whole basis for phase 2).
4. **The config's linked docs** — fetch `links.documentation` (settings / API pages give the
   real parameter list, defaults, and quality-vs-speed notes) and `links.repository`.
5. **The actual publication — verify, do not trust `references.doi` blindly.** The DOI in the
   config is often the *generic/framework* paper, not the one describing the model actually
   run (real example: `cellposev4`'s DOI pointed at the original 2021 Cellpose paper, not the
   Cellpose-SAM preprint it runs). Find and READ the method-specific paper; cross-check claims
   against its abstract rather than paraphrasing a README or search snippet.
   - bioRxiv/medRxiv **block WebFetch (403)** → use the bioRxiv MCP `get_preprint` tool by DOI
     for title/authors/abstract/metadata. For journals, WebFetch the DOI (follow the redirect
     to the publisher host on the second call) or use the PubMed MCP.
   - If the config DOI is wrong or generic, **flag it and offer to fix** — cite both the
     method-specific paper *and* the framework paper (Viash `references.doi` accepts a list).

**Write it down in `NOTES.md`.** Use the **`document-method-troubleshooting`** skill for the
NOTES.md + memory mechanics and section conventions — do not reinvent them here. Add a
**"Citation caveat"** near the top if the DOI didn't match, and record which defaults are
speed-tuned (you'll need that in phase 2).

## Phase 2 — Identify optimization parameters + criteria

Enumerate **every** knob the underlying tool's `eval`/API exposes (from phase 1's docs), then
organise them. The organising principle is *impact tiers* — this IS the criteria:

- **Tier 0 — the input, not a parameter.** Often the biggest lever: channels fed, resolution,
  which stain. (e.g. feeding a 2-channel membrane+nuclear stack vs a single channel.) Call it
  out even though it isn't in the sweep.
- **Tier 1 — highest impact on output quality.** The knobs that change *what* gets detected
  (recall↔precision dials, size/scale params).
- **Tier 2 — quality/speed trade-offs.** Knobs that mostly buy accuracy at a time cost.
- **Tier 3 — not exposed by the component.** High-value knobs worth *adding* to the config
  for a serious sweep.

Three rules that come straight from phase 1:

- **Defaults are usually speed-tuned**, so "optimize for quality" mostly means *walking them
  back toward the tool's own defaults* (re-enable a QC check that was disabled for speed, raise
  an iteration count that was floored, etc.). Say so.
- **Ground every default and range in evidence** (docs + paper), not memory. Tie ranges to the
  data when you can (e.g. pick a `diameter` range from the pixel size: Xenium ~0.2125 µm/px ⇒
  nuclei ~35–40 px, whole cells ~60–70 px).
- **Any parameter already set to a non-default value is a sweep candidate, not a silent
  baseline.** Compare each config default (and anything the run script's `default:` block sets)
  against the *tool's* true default. A deviation is evidence someone already found the knob
  matters — so it must go on a sweep axis, with its range straddling both the deviated value
  and the tool default. The **only** exception is strictly-performance / resource knobs that
  can't change the output (e.g. cores/threads, `batch_size`, `n_workers`, memory, tile size for
  pure memory management) — those stay fixed in `default:` and never enter the sweep.

**Criteria for a good sweep axis** (used in phase 3):

- Continuous / multi-level params get a **range of several values**.
- **Boolean** params have exactly **one** meaningful non-default value (the flip) — don't pad
  them.
- **Omit the default value from the sweep list** — the "default variant" already covers that
  point (see phase 3); repeating it just runs a duplicate.

**If a Tier-1/Tier-3 knob isn't exposed, expose it.** Add it to `config.vsh.yaml`'s
`arguments:` (type, default = the tool's default, a clear description) and thread it into the
script's eval-params (match the existing pattern — e.g. add the key to the params-forwarding
comprehension). Validate with `viash config view <config>`. Record the new arg and the tier
ranking in NOTES.md under an **"Optimization / tuning"** section.

## Phase 3 — Set up the sweep run scripts

### How the sweep expands (know this cold)

The `run_benchmark` workflow reads a `method_parameters_yaml` file
(`src/workflows/run_benchmark/main.nf`, init map ~L22-32) shaped as:

```yaml
parameters:
  <component>:
    default: { argA: ..., argB: ... }   # baseline args applied to EVERY variant
    sweep:   { argA: [v1, v2, v3] }      # one extra variant per value, that ONE arg overridden
```

The expansion (main.nf ~L738-773) is a **star around the default, NOT a grid**: sweep list
lengths **add**, they do not multiply. **Total variants = 1 (default) + Σ(len of each sweep
list).** There's also a hard constraint — **at most one non-default method/variant at a time**
across the whole pipeline — so the swept method must be the only non-default thing; keep every
other stage on its single default.

### Build the script

Cloud/nebius only — **do NOT generate a local `run_test_<method>_local.sh`**. Everything the
sweep needs lives under **`scripts/run_benchmark/param_sweep/`** (create the folder if missing):

- **`param_sweep/<method>_params.yaml`** — the `method_params` (default + sweep) as a committed
  file, NOT an inline heredoc. This is the single source of truth for the sweep.
- **`param_sweep/run_test_<method>_nebius.sh`** — copy `run_test_nebius.sh` / `run_gpu_nebius.sh`
  (`tw launch`); GPU methods also keep the `gpu` label + GPU compute env.

In the settings block: enable **all default methods + the target method** at its stage
(comment the rest), and set `method_parameters_yaml` to the committed params file's raw-GitHub
URL (see below).

### The params-file placement gotcha (cloud)

`readYaml` uses Nextflow's `file()` (`common/nextflow_helpers/workflowHelper.nf`), and the
workflow opens `method_parameters_yaml` **at runtime on the cloud compute env**, so a local
path (e.g. `/tmp`) does not exist there — this is why the binning `method_params` block is
commented out in the stock `run_test_nebius.sh`. Note `/scratch` is **read-only from the launch
host** (the cloud writes results there, you can't stage into it), so it is NOT a staging
option. `file()` does stage **`http(s)://`** and **`s3://`**, so use one of:

- **Committed file + raw GitHub URL (preferred).** `tw launch` already pulls the repo from
  GitHub, so commit `scripts/run_benchmark/param_sweep/<method>_params.yaml` and set
  `method_parameters_yaml` to its
  `https://raw.githubusercontent.com/<org>/<repo>/<branch>/scripts/run_benchmark/param_sweep/<method>_params.yaml`
  URL (public repo ⇒ no auth). The file must be **pushed to `<branch>` before launching**;
  `curl -fsSL` it first in the script to fail fast if it's not there yet. This is independent
  of `--revision` (which selects the pipeline *code*).
- **`s3://` URI** — upload the params to a bucket the launcher can write and the compute env
  can read (`aws s3 cp`, possibly `--profile`).

GPU methods also need the `gpu` label + GPU compute env.

### Validate before handing back

- `bash -n param_sweep/run_test_<method>_nebius.sh` (syntax).
- **Non-default audit (do this before finalizing the sweep).** Walk every argument in
  `config.vsh.yaml` and every key in the params `default:` block, and compare its effective
  value to the *tool's* true default (from phase 1). For each one that deviates:
  - if it's a strictly-performance / resource knob (cores, threads, `batch_size`, `n_workers`,
    memory, memory-only tiling) → leave it fixed and note it;
  - otherwise → it **must** appear as a `sweep:` axis. If it doesn't, add it (range straddling
    the deviated value and the tool default, per the phase-2 rule) before handing back.
  List what you found and how you resolved it.
- Run `validate_sweep.py` (in this skill dir) on `param_sweep/<method>_params.yaml` to print the
  per-param value counts, the **total variant count**, and a warning if any sweep list repeats
  its default or if the sweep is heavy for a "test" run.
- `curl -fsSL` the raw-GitHub params URL to confirm it resolves once the file is pushed.
- Sanity-check the variant count against the run's budget — a 19-variant sweep × the full
  downstream pipeline is a lot for a quick smoke test; note it and offer to trim.

## Output

Tell the user, briefly: the NOTES.md sections written, any config args newly exposed (and that
they need `viash ns build` + a container rebuild to take effect — see `check-component`), the
`param_sweep/` files created (params yaml + nebius script), and the resulting variant count.
Don't reprint whole files.
