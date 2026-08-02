# split — NOTES

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the `split`
expression-correction method. Committed with the code.

## What this component is

- **Stage / API:** `methods_expression_correction` (merges
  `/src/api/comp_method_expression_correction.yaml`). Runs late in the iST pipeline: it takes
  `spatial_with_cell_types.h5ad` (cell x gene matrix with cell-type labels, centroids, a
  `counts` and a `normalized` layer) plus the `scrnaseq_reference.h5ad`, and rewrites the
  spatial object's `counts` / `normalized` layers with purified (contamination-corrected) ones.
  Unlike resolvi_correction, split **requires the SC reference** (it deconvolves against it).
- **Underlying tools:** **SPLIT** (`bdsc-tds/SPLIT`, R) layered on top of **RCTD** (spacexr).
  RCTD deconvolves each spatial cell against the SC reference in `doublet_mode`; SPLIT then
  "purifies" the counts — for a doublet it splits the mixed transcripts to the primary cell
  type, removing spillover/contamination. R method (`script.R`), CPU-only.
- **GPU:** none. Nextflow directives `[ hightime, highcpu, highmem ]`; no `gpu` label.
- **Links:** docs/repo https://github.com/bdsc-tds/SPLIT .

## Citation

`references.doi: 10.1101/2025.04.23.649965` is **correct** (verified via bioRxiv): Bilous,
Buszta, et al., 2025, *"From Transcripts to Cells: Dissecting Sensitivity, Signal
Contamination, and Specificity in Xenium Spatial Transcriptomics"* — the paper that introduces
SPLIT (Spatial Purification of Layered Intracellular Transcripts). Method-specific, not a
generic framework citation. No caveat needed.

## script.R — step by step

1. Read `input_spatial_with_cell_types` twice — as `SingleCellExperiment` (`sce`) and as a
   `Seurat` object (`xe`) (L30-31). `sce` carries the assays/coords; `xe` supplies the raw
   `counts` layer handed to `SPLIT::purify`.
2. If `!keep_all_cells`: drop zero-count cells from both objects (L34-38).
3. Build the `SpatialRNA` puck from `centroid_x/centroid_y` + the `counts` assay (L41-46).
4. Read the SC reference; **drop reference cells with zero panel-gene counts** (L56-60) — after
   subsetting to the shared ~few-hundred-gene panel, some cells have `nUMI==0`, which makes
   spacexr's per-cell normalisation `NaN` and poisons every cell type's DE-gene means (see the
   memory note `rctd-split-zero-umi-reference-nan`). This guard is load-bearing.
5. Keep only cell types with >=25 cells (RCTD minimum) (L62-64); sanitise `/` out of cell-type
   names (spacexr's `check_cell_types` rejects them) (L76); build the `Reference(min_UMI=0)`.
6. `create.RCTD(puck, reference, max_cores, gene_cutoff, fc_cutoff, gene_cutoff_reg,
   fc_cutoff_reg, UMI_min, UMI_min_sigma)` (L89-94) — **this is where all six exposed threshold
   args are forwarded.** Then `run.RCTD(doublet_mode = "doublet")` (L95) — `doublet_mode` is
   **hard-coded**.
7. `SPLIT::run_post_process_RCTD(myRCTD)` (L106) then `SPLIT::purify(counts=<Seurat counts>,
   rctd=RCTD, DO_purify_singlets=TRUE)` (L110-114) — `DO_purify_singlets` is **hard-coded TRUE**
   and `DO_remove_residual_contamination` is **left at its package default (not set)**.
8. Stash the incoming `normalized` as `normalized_uncorrected` (L121); write purified counts
   into `corrected_counts` (copy counts, then overwrite the updated cells) (L124-127);
   `logNormCounts` with library-size factors -> new `normalized` (L130-131); `write_h5ad` (L136).

## Arguments (exposed in config.vsh.yaml)

All six thresholds map straight into the `create.RCTD(...)` call; `keep_all_cells` is a
component-local zero-cell guard (no upstream equivalent).

| arg | type | config default | spacexr (RCTD) default | notes |
|-----|------|----------------|------------------------|-------|
| `--keep_all_cells` | boolean | `false` | n/a | drop zero-count cells; description warns TRUE "may cause errors". Robustness flag, not a tuning lever. |
| `--gene_cutoff` | double | `0.0` | **0.000125** | min normalized mean expr, platform-effect DE gene |
| `--fc_cutoff` | double | `0.1` | **0.5** | min log-FC, platform-effect DE gene |
| `--gene_cutoff_reg` | double | `0.0` | **0.0002** | min normalized mean expr, regression DE gene |
| `--fc_cutoff_reg` | double | `0.1` | **0.75** | min log-FC, regression DE gene |
| `--umi_min` | integer | `20` | **100** | min total UMI per spatial cell to include |
| `--umi_min_sigma` | integer | `20` | **300** | min UMI for cells fitting the platform-effect variance |

**Why the six defaults are relaxed (deliberate, not accidental):** RCTD's DE-gene / UMI
thresholds are tuned for whole-transcriptome references (~20k genes). iST panels are small
(~100-500 genes) and low-count, so the stock spacexr thresholds strand too few DE genes and
`create.RCTD` aborts with *"fewer than 10 regression differentially expressed genes"*. These
defaults mirror the standalone `rctd` cell-type-annotation component exactly.

**Hard-coded (NOT exposed) SPLIT/RCTD knobs that matter:** `doublet_mode="doublet"` (L95),
`DO_purify_singlets=TRUE` (L113), `DO_remove_residual_contamination` (unset -> package default,
L110-114), `Reference(min_UMI=0)` (L78), the `>=25`-cells-per-type filter (L63).

## Optimization / tuning

**Verified upstream defaults** — spacexr `create.RCTD`: `gene_cutoff=0.000125`, `fc_cutoff=0.5`,
`gene_cutoff_reg=0.0002`, `fc_cutoff_reg=0.75`, `UMI_min=100`, `UMI_min_sigma=300`. SPLIT
`purify` (from the repo README): `DO_purify_singlets` and `DO_remove_residual_contamination`
(the latter "removes an additional 2-5% of counts" for improved specificity, v0.3.0+).

**Non-default audit result:** all six exposed thresholds sit at RELAXED, non-default values
(walked away from the spacexr defaults). None is a performance/resource knob, so per the
phase-3 audit **each is forced onto a sweep axis** — the deviation is itself evidence the knob
matters. `keep_all_cells` is the one exposed arg NOT swept: it has no upstream tool default to
straddle and its own description warns that flipping it to TRUE "may cause errors", so it stays
fixed at `false`. (This is the "do not sweep every knob" — 6 of the 7 exposed args are swept.)

Tiers:

- **Tier 0 — input, not a parameter.** The SC reference quality/annotation and the shared panel
  gene set (subsetting to a small panel is what forces the relaxed thresholds in the first
  place); upstream segmentation/assignment quality. Biggest levers overall.
- **Tier 1 — the exposed thresholds (what the sweep varies).** Each relaxed default is walked
  back TOWARD its spacexr default (the meaningful direction; the relaxed default already sits at
  the permissive extreme):
  - `fc_cutoff`: 0.1 -> **0.25, 0.5**   (sharpest — platform-effect DE-gene selection)
  - `fc_cutoff_reg`: 0.1 -> **0.4, 0.75**   (sharpest — regression DE-gene selection)
  - `gene_cutoff`: 0.0 -> **0.0000625, 0.000125**   (companion expression-mean filter)
  - `gene_cutoff_reg`: 0.0 -> **0.0001, 0.0002**   (companion expression-mean filter)
  - `umi_min`: 20 -> **50, 100**   (which spatial cells get deconvolved)
  - `umi_min_sigma`: 20 -> **100, 300**   (which cells fit the platform-effect variance)
  Total variants = 1 default + 6 axes x 2 values = **13**.
- **Tier 3 — high-value SPLIT-specific knobs NOT surfaced by the component (future work; would
  need a config arg + `viash ns build` + container rebuild — NOT in the submittable sweep):**
  - `DO_purify_singlets` (SPLIT `purify`, hard-coded TRUE, script.R L113): whether to also
    purify singlet cells, not just doublets. Arguably the #1 SPLIT-specific lever. Expose as a
    boolean (`--do_purify_singlets`, default true) and sweep the flip to `false`.
  - `DO_remove_residual_contamination` (SPLIT `purify`, unset, v0.3.0+): removes an extra 2-5%
    of non-reference-supported counts for higher specificity. Expose as a boolean
    (`--remove_residual_contamination`, default false) and sweep the flip to `true`.
  - `doublet_mode` (RCTD `run.RCTD`, hard-coded "doublet", script.R L95): "doublet" is what
    SPLIT's purification consumes; "full"/"multi" change behaviour materially — expose only if
    the purification path is adapted to accept them.

  Exposing any of these needs: add to `config.vsh.yaml arguments:`, thread `par[...]` into the
  matching `SPLIT::purify(...)` / `run.RCTD(...)` call in `script.R`, then `viash ns build` + a
  container rebuild (see `check-component`) before a build/main Nebius run would honour them.

## Wiring

- Registered in the run_benchmark workflow: dep list
  `src/workflows/run_benchmark/config.vsh.yaml` L186; present in the `--expression_correction_methods`
  default string L107 (`no_correction:resolvi_correction:split`).
- Sweep: `scripts/run_benchmark/param_sweep/split_params.yaml` (committed params) +
  `scripts/run_benchmark/param_sweep/run_test_split_nebius.sh` (Nebius launcher, standard CPU
  compute env `5hfmdCBxMRd4nHZaJKYEQZ`, labels `task_ist_preprocessing,test,split` — no `gpu`).
  The launch keeps `no_correction` enabled alongside `split` so the sweep is scored against the
  uncorrected baseline; every other stage stays on its single default.

## Risk points / gotchas

- **Submittable-now constraint:** the Nebius run uses `--revision build/main`, so the sweep can
  only vary args already in the build/main container (the six thresholds). Any newly-exposed
  Tier-3 arg is silently ignored until rebuilt.
- **Zero-UMI reference-cell guard (script.R L56-60) is load-bearing** — without it small-panel
  references abort RCTD with "fewer than 10 DE genes" (memory: `rctd-split-zero-umi-reference-nan`).
- **High-threshold sweep variants can legitimately fail** — pushing a threshold back to the
  spacexr default on a very small panel can drop under the 10-DE-gene floor and abort inside
  `create.RCTD`. That failure boundary is part of what the sweep measures, not a bug.
- **SPLIT installed from HEAD** (`remotes::install_github("bdsc-tds/SPLIT")`, unpinned) — the
  `purify(counts=, rctd=, ...)` convenience signature the script relies on tracks whatever is on
  the default branch; a breaking upstream API change would surface only at container rebuild.
