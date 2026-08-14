# pciseq (transcript assignment) — NOTES

`src/methods_transcript_assignment/pciseq/` implements
`src/api/comp_method_transcript_assignment.yaml`. It wraps **pciSeq** (Qian et al. 2020,
*Nat Methods*, DOI `10.1038/s41592-019-0631-4` — the correct pciSeq paper; no citation
caveat), installed from `acycliq/pciSeq@issue_14`, via `txsim.preprocessing.run_pciSeq`.

pciSeq is a **Variational-Bayes probabilistic model** that jointly (a) cell-types each
segmented cell against an scRNA-seq reference and (b) assigns every transcript ("spot") to a
candidate parent cell or to a background/misread class. It therefore needs three inputs: the
transcripts (points), a **segmentation label image** (cells come from an upstream segmenter —
here `custom_segmentation`), and an **scRNA-seq reference** (`input_scrnaseq`, counts layer).
CPU-only.

## How the args reach the tool

`script.py` L176-188 forwards nine `opts` keys straight into pciSeq (`opts_keys` ->
`opts = {k: par[k] ...}`), and `txsim.run_pciSeq` passes that dict unmodified to
`pciSeq.fit(..., opts=opts)` (it only auto-fills `exclude_genes` for panel genes missing from
the reference). So all nine exposed knobs are genuinely live:
`max_iter, CellCallTolerance, rGene, Inefficiency, InsideCellBonus, MisreadDensity, SpotReg,
nNeighbors, rSpot, save_data`. `dtype` is present in the VIASH-START stub but is NOT in
`opts_keys`, so it is not forwarded (memory-only knob anyway). `transcripts_key` /
`coordinate_system` are structural. Several load-bearing fixes live in script.py (assertion
monkey-patch for pciSeq issue #14, duplicate-index reset before `sd.transform`, out-of-bounds
spot pre-filter on both borders, all-zero-segmentation guard — see the inline comments and the
`pciseq-empty-segmentation-lung-merscope` memory).

## Optimization / tuning

**NON-DEFAULT AUDIT (result: clean).** Every one of the nine exposed `opts` defaults in
`config.vsh.yaml` equals pciSeq's own `config.py` default **exactly** (verified against
`acycliq/pciSeq@issue_14`): `max_iter=1000`, `CellCallTolerance=0.02`, `rGene=20`,
`Inefficiency=0.2`, `InsideCellBonus=2`, `MisreadDensity=1e-5`, `SpotReg=0.1`, `nNeighbors=3`,
`rSpot=2`, `save_data=False`. Nothing was walked away from the tool default, so **no knob is
forced onto a sweep axis by the "deviation is evidence" rule** — ranges are chosen on first
principles and each straddles pciSeq's default. (This is unlike the sibling RCTD/SPLIT sweeps,
whose defaults are all deliberately relaxed.)

Tiers:

- **Tier 0 — input, not a parameter (biggest levers overall).** The scRNA-seq reference
  quality/annotation, the shared gene panel, and — because pciSeq assigns spots *relative to
  cell boundaries* — the upstream **segmentation** it is handed. A poor label image caps what
  any `opts` tuning can recover.
- **Tier 1 — what gets assigned to which cell (swept).** Each straddles pciSeq's default:
  - `InsideCellBonus`: 2 -> **0, 1, 4, 6**  (sharpest — log-likelihood boost for a spot inside a
    cell boundary; high = trust segmentation more, 0/low = let expression reassign spots)
  - `MisreadDensity`: 1e-5 -> **1e-6, 1e-4, 1e-3**  (background/misread rate per pixel; high =
    more spots called background => sparser/higher-precision, low = assign nearly everything)
  - `nNeighbors`: 3 -> **2, 6**  (candidate parent cells per spot; a background super-neighbor is
    always added on top)
  - `Inefficiency`: 0.2 -> **0.1, 0.5**  (mean of the per-gene efficiency Gamma prior:
    observed/theoretical count ratio)
  - `rGene`: 20 -> **10, 40**  (variance of that efficiency Gamma prior; higher = tighter)
  Total variants = 1 default + (4 + 3 + 2 + 2 + 2) = **14**.
- **Tier 2 — exposed but deferred from this first (light) test run.** All live knobs, add to the
  sweep for a deeper pass: `SpotReg` (additive expression-mean regularizer, default 0.1),
  `rSpot` (Negative-Binomial count dispersion, default 2), `CellCallTolerance` (0.02) and
  `max_iter` (1000) — the last two are a convergence budget (quality/speed), least likely to move
  final quality, so they are fixed here.
- **Tier 3 — pciSeq priors NOT surfaced by the component (future work; need a new `config.vsh.yaml`
  arg + thread into `opts_keys` + `viash ns build` + container rebuild).** From pciSeq
  `config.py`: `cell_radius` (None), `cell_type_prior` ('uniform'),
  `mean_gene_counts_per_class` (60), `mean_gene_counts_per_cell` (15). Not exposed because the
  sharpest lever (`InsideCellBonus`) already is; these are lower-value priors.

**SUBMITTABLE AS-IS.** No config/script edits were made — all five swept args already exist in
`config.vsh.yaml`, so the `--revision build/main` Nebius run honours them with no rebuild.

## Wiring

- Registered in the run_benchmark workflow: dep list
  `src/workflows/run_benchmark/config.vsh.yaml` L166; present in the
  `--transcript_assignment_methods` default string L71
  (`basic_transcript_assignment:baysor:clustermap:pciseq:comseg:proseg:segger:fastreseg`).
- Sweep: `scripts/run_benchmark/param_sweep/pciseq_params.yaml` (committed params) +
  `scripts/run_benchmark/param_sweep/run_test_pciseq_nebius.sh` (Nebius launcher, standard CPU
  compute env `5hfmdCBxMRd4nHZaJKYEQZ`, labels `task_ist_preprocessing,test,pciseq` — no `gpu`).
  The launch keeps `basic_transcript_assignment` enabled alongside `pciseq` so the sweep is
  scored against the baseline assignment; every other stage stays on its single default.

## Risk points / gotchas

- **Params-file placement (cloud):** `method_parameters_yaml` is opened by the WORKFLOW at
  runtime on the compute env, so a local `/tmp` path won't exist and `/scratch` is read-only from
  the launch host. The sweep is therefore a COMMITTED file read from GitHub via a raw URL —
  commit AND push `pciseq_params.yaml` to `$params_branch` before launching (independent of
  `--revision`, which selects pipeline CODE).
- **Submittable-now constraint:** the run uses `--revision build/main`, so any newly-exposed
  Tier-3 arg is silently ignored until `viash ns build` + a container rebuild land in build/main
  (see `check-component`).
- **Empty / all-zero segmentation** aborts pciSeq with a cryptic "zero-size array to reduction
  maximum"; script.py raises an actionable error instead (memory:
  `pciseq-empty-segmentation-lung-merscope`). High-`InsideCellBonus` / low-`MisreadDensity`
  variants push more spots into cells but cannot conjure cells that segmentation didn't produce.
