# RCTD — cell type annotation (NOTES)

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the RCTD
component. The running iteration log lives in the memory file
`rctd-split-zero-umi-reference-nan.md` (repo auto-memory) — that pointer + log, this
file the depth.

## What this component is

- **Stage / API:** cell type annotation (`src/api/comp_method_cell_type_annotation.yaml`).
  Inputs `input_spatial_normalized_counts` (the aggregated cell x gene `.h5ad` with
  centroids) and `input_scrnaseq_reference` (labelled scRNA-seq atlas); outputs the same
  spatial object with a `cell_type` column in `.obs`/`colData`.
- **Unusual:** this is an **R** component (`script.R`, `openproblems/base_r:1`) that wraps
  **spacexr::RCTD** (Robust Cell Type Decomposition). It installs `spacexr` from GitHub
  (`dmcable/spacexr`) at build time — an un-pinned `install_github`, so the exact spacexr
  version floats with the build date (last observed live: 2.2.1).
- **Links:** docs/repo `https://github.com/dmcable/spacexr`.
- **Publication:** Cable et al., *Robust decomposition of cell type mixtures in spatial
  transcriptomics*, Nat Biotechnol 2022 (DOI `10.1038/s41587-021-00830-w`, PMID 33603203).
  **The config DOI is correct** (it is the RCTD method paper, not a generic framework
  paper) — no citation caveat. The paper's headline contribution is the `doublet_mode`
  decomposition (each pixel assigned at most two cell types), and platform-effect / DE-gene
  selection is the machinery the exposed thresholds control.

## script.R — step by step

1. **L27** read the spatial `.h5ad` as a `SingleCellExperiment`.
2. **L30-35** build a `SpatialRNA` "puck" from `centroid_x/centroid_y` and the raw
   `counts` assay (RCTD works on **raw integer counts**, not the normalized layer — despite
   the input being named `*_normalized_counts`, the raw `counts` assay is what is read).
3. **L38** read the scRNA-seq reference as a `SingleCellExperiment`.
4. **L45-49 — zero-UMI reference fix (load-bearing).** After `process_dataset` subsets the
   reference to the shared spatial panel (~hundreds of genes), some reference cells express
   none of those genes (`nUMI == 0`). spacexr's `get_cell_type_info` divides each cell by
   its `nUMI` → `NaN`; a single `NaN` column poisons `other_mean` (rowMeans) for **every**
   cell type in `get_de_genes` → all logFC `NaN` → 0 DE genes → `create.RCTD` aborts with
   "fewer than 10 regression differentially expressed genes". Dropping zero-count cells here
   is the fix (see memory `rctd-split-zero-umi-reference-nan.md`). **This — not the
   thresholds — was the true root cause of the historical abort.**
5. **L52-53 — CELL_MIN_INSTANCE, applied by hand.** Keep only reference cell types with
   >=25 cells. This replicates RCTD's `CELL_MIN_INSTANCE=25` default but is a **hardcoded
   pre-filter, not an exposed arg** (so it never reaches `create.RCTD`).
6. **L82-87** sanitize cell-type factor levels containing `/` (spacexr's
   `check_cell_types` rejects them); keep a safe→original name map to restore labels later.
7. **L89** `Reference(ref_counts, cell_types, min_UMI = 0)` — `min_UMI` hardcoded to 0 so
   the (now zero-UMI-free) cells all pass.
8. **L100-105** `create.RCTD(...)` — the **6 exposed threshold args** are forwarded here:
   `gene_cutoff`, `fc_cutoff`, `gene_cutoff_reg`, `fc_cutoff_reg`, `UMI_min`,
   `UMI_min_sigma`. `max_cores` = `meta$cpus` (performance only).
9. **L106** `run.RCTD(myRCTD, doublet_mode = "doublet")` — `doublet_mode` is **hardcoded**,
   not exposed.
10. **L109-114** take `results_df$first_type`; pixels with `spot_class == "reject"` are
    relabelled `None_sp`.
11. **L117-128** write predictions back into `colData(sce)$cell_type` (restoring the
    original `/`-containing names) and `write_h5ad` the result.

## Arguments

Six exposed args, all `create.RCTD` DE-gene / UMI thresholds. **Every config default is
deliberately relaxed away from the spacexr default** — RCTD's defaults are tuned for
whole-transcriptome references (~20k genes) with high per-cell UMIs, whereas iST uses small
curated panels (~100-500 genes) with compressed fold-changes and low per-cell counts, so the
stock defaults strand too few DE genes / drop too many cells.

| Arg (config) | Config default | spacexr default | Maps to | Meaning |
|---|---|---|---|---|
| `--gene_cutoff` | 0.0 | 0.000125 | `create.RCTD(gene_cutoff)` | min normalized mean expr for a platform-effect DE gene |
| `--fc_cutoff` | 0.1 | 0.5 | `create.RCTD(fc_cutoff)` | min log-FC for a platform-effect DE gene |
| `--gene_cutoff_reg` | 0.0 | 0.0002 | `create.RCTD(gene_cutoff_reg)` | min normalized mean expr for a regression DE gene |
| `--fc_cutoff_reg` | 0.1 | 0.75 | `create.RCTD(fc_cutoff_reg)` | min log-FC for a regression DE gene |
| `--umi_min` | 20 | 100 | `create.RCTD(UMI_min)` | min total UMI per spatial cell to annotate it |
| `--umi_min_sigma` | 20 | 300 | `create.RCTD(UMI_min_sigma)` | min UMI for cells used to fit platform-effect variance |

Hardcoded (NOT exposed): `doublet_mode="doublet"` (L106), the `>=25` cells/type filter
(≈`CELL_MIN_INSTANCE`, L52), `Reference(min_UMI=0)` (L89).

## Setup / Docker

Base `openproblems/base_r:1`; installs `SingleCellExperiment`, `anndataR`, `rhdf5`,
`devtools` via Bioc, then `devtools::install_github('dmcable/spacexr', build_vignettes=FALSE)`.
The `SingleCellExperiment` reinstall comment (config L63-66) is a workaround for a
`SpatialExperiment`/Seurat install-order bug. **The `install_github` is not commit-pinned**,
so the spacexr version floats — a risk point if upstream changes DE-gene defaults or the
`create.RCTD`/`run.RCTD` signatures.

## Wiring

- Registered as `rctd` in the `celltype_annotation_methods` stage of the run_benchmark
  workflow. Stage default is `tacco`.
- Nextflow labels (config L75): `hightime, midcpu, highmem`.
- Sweep scripts: `scripts/run_benchmark/param_sweep/rctd_params.yaml` +
  `run_test_rctd_nebius.sh` (see below).

## Risk points / gotchas

- **Small-panel DE-gene floor.** Even with the zero-UMI fix, raising the fold-change
  thresholds toward the spacexr defaults shrinks the DE-gene set; on a very small panel this
  can drop back under the 10-gene floor and re-trigger the `create.RCTD` abort. The relaxed
  defaults exist precisely to stay clear of that floor. **Consequence for the sweep: the
  high-threshold variants (values approaching the spacexr defaults) may legitimately fail on
  small-panel datasets** — that failure boundary is part of what the sweep measures.
- Un-pinned spacexr (version floats with build date).
- `doublet_mode` and the `CELL_MIN_INSTANCE`-equivalent filter are hardcoded — see
  Optimization / tuning below.

## Optimization / tuning

Impact tiers (grounding: Cable et al. 2022 + spacexr `create.RCTD`/`run.RCTD` docs). The
sweep in `param_sweep/rctd_params.yaml` varies only **already-exposed** args (so it is
submittable against the current build/main container without a rebuild).

- **Tier 0 — input, not a parameter.** The scRNA-seq **reference**: its cell-type
  granularity and, above all, its **overlap with the spatial gene panel**. RCTD learns
  profiles only on the shared-panel genes, so panel size/quality dominates everything below.
  Not a sweep axis.
- **Tier 1 — highest impact on quality (EXPOSED, all six on the sweep).** The DE-gene /
  UMI thresholds. `fc_cutoff` / `fc_cutoff_reg` are the sharpest levers (they set which genes
  are "informative"); `umi_min` sets which cells get annotated at all vs dropped;
  `gene_cutoff` / `gene_cutoff_reg` / `umi_min_sigma` are companion filters. **All six are
  set to a non-default (relaxed) value, so per the non-default rule all six are on a sweep
  axis, each range straddling the relaxed config default and the spacexr default.** Because
  the relaxed default already sits at/near the permissive extreme, the meaningful sweep
  direction is *toward* the stricter spacexr default.
- **Tier 1 but NOT exposed → Tier 3 (deferred).** `doublet_mode` (`run.RCTD`, hardcoded
  `"doublet"`). This is RCTD's single biggest behavioural lever: `"doublet"` assigns ≤2
  types/pixel (the paper's headline mode, aimed at mixed Slide-seq pixels), `"full"` does
  unrestricted deconvolution, `"multi"` is a greedy multi-type extension. For one-cell-per-
  object segmented iST, `doublet` + `first_type` is a reasonable default, but `full`/`multi`
  could change results materially. **Worth exposing** as `--doublet_mode` (thread into the
  `run.RCTD` call). Not in this sweep — exposing it needs `viash ns build` + a container
  rebuild, so a stale build/main container would silently ignore it.
- **Tier 3 (deferred) — other un-exposed knobs.** `CELL_MIN_INSTANCE` (currently the
  hardcoded `>=25` cells/type pre-filter, L52) — exposing it would let the sweep trade rare-
  cell-type coverage against profile stability. `Reference(min_UMI=0)` (L89) is intentional
  given the zero-UMI fix and is best left fixed.
- **Performance only (fixed, never swept):** `max_cores` (= `meta$cpus`).

**To promote any Tier-3 knob into a real sweep:** add the arg to `config.vsh.yaml`
(`type`, `default` = the spacexr default), thread it into the `run.RCTD`/`create.RCTD` call
in `script.R`, `viash config view` to validate, then `viash ns build` + rebuild the
container (see `check-component`) before it can be swept on the cloud.
