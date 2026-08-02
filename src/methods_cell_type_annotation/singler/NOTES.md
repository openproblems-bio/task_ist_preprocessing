# singler — developer notes

## What this component is

Cell-type annotation method (stage API `src/api/comp_method_cell_type_annotation.yaml`,
subtype `method_cell_type_annotation`). It is a **thin Python wrapper around
[singler-py](https://github.com/SingleR-inc/singler-py)** (the Python port of the
Bioconductor `SingleR`), a **reference-correlation labeller**: it builds per-label marker
sets from an annotated scRNA-seq reference, then assigns each spatial cell the label whose
reference profile its expression correlates with best (Spearman, per-label score at a
correlation quantile), optionally refined by fine-tuning.

- CPU-only (no GPU); resource label `[ midtime, midcpu, midmem ]`.
- Docker: `openproblems/base_python:1` + `pypi: [singler]`, merging
  `/src/base/setup_spatialdata_partial.yaml`. Merges the base setup and needs no version
  gymnastics — nothing load-bearing to record here.
- Links: docs/repo `https://github.com/SingleR-inc/singler-py`.

### Citation note
`references.doi: 10.1038/s41590-018-0276-y` is **correct** — Aran et al., *Nat Immunol*
2019, the paper that introduced SingleR. (Unlike some components in this repo, the DOI is
the method-specific paper, not a generic framework paper, so no fix is needed.)

## script.py — step by step

1. `sce.read_h5ad(input_spatial_normalized_counts)` and a parallel `ad.read_h5ad` of the
   same file (L27-28). The `SingleCellExperiment` is used for the matrix; the `AnnData` is
   the object we write labels back onto.
2. `sce_ref = sce.read_h5ad(input_scrnaseq_reference)` (L30).
3. Test matrix = the spatial **`counts`** assay, `mat.sorted_indices()` (L34-35). Reference
   matrix = the reference **`normalized`** assay, `sorted_indices()` (L37-38).
   `sorted_indices()` is the "magic line" that puts the CSR/CSC into the layout singler
   expects.
4. `singler.train_single(ref_data=mat_ref, ref_labels=<ref cell_type column>,
   ref_features=..., test_features=...)` (L41-44) — builds the prebuilt reference.
5. `singler.classify_single(mat, ref_prebuilt=built)` (L47) — classifies.
6. `adata_sp.obs["cell_type"] = output["best"]` (L49), then `adata_sp.write(output)` (L53).

### Two things that bite (read before touching this)

- **The two config-exposed non-IO args are DEAD.** `--celltype_key` (default `cell_type`)
  and `--labels_key` (default `cell_labels`) are declared in the merged config but the
  script **never reads `par['celltype_key']` or `par['labels_key']`**. `ref_labels` is
  hardcoded to `sce_ref.get_column_data().column("cell_type")` (L42), and the output column
  is hardcoded to `"cell_type"` (L49). So changing either arg on the command line (or via a
  sweep) currently has **no effect whatsoever**. See "Optimization / tuning" — wiring
  `celltype_key` is the single highest-leverage fix and the basis of the prepared sweep.
- **Test matrix is raw `counts`, reference is `normalized`.** SingleR scores with a
  rank-based (Spearman) correlation, which is invariant to per-cell monotonic scaling, so a
  raw-counts test side is mostly harmless; but marker detection on the reference does expect
  log-normalized input, which is why the reference uses `normalized`. Not a tuning axis,
  just context.

## Arguments

| arg | default | what it maps to upstream | status |
|-----|---------|--------------------------|--------|
| `--input_spatial_normalized_counts` | — | test matrix (`counts` assay) + AnnData to annotate | used |
| `--input_scrnaseq_reference` | — | `train_single(ref_data=…)` (`normalized` assay) | used |
| `--input_transcript_assignments` | — | (optional, unused by script) | ignored |
| `--celltype_key` | `cell_type` | *should* select the reference label column for `ref_labels` | **declared but NOT wired** (hardcoded `"cell_type"`) |
| `--labels_key` | `cell_labels` | (spatial label key) | **declared but NOT wired** (never read) |
| `--output` | — | `adata_sp.write(...)` | used |

None of singler-py's algorithm knobs (`marker_method`, `num_de`, `quantile`,
`use_fine_tune`, `fine_tune_threshold`, `aggregate`) are exposed — `train_single` /
`classify_single` are called entirely at their upstream defaults.

## Wiring

- Registered in `src/workflows/run_benchmark/config.vsh.yaml`: dependency list
  (`methods_cell_type_annotation/singler`) and the `--celltype_annotation_methods` default
  string `ssam:tacco:moscot:mapmycells:tangram:singler:rctd`; fan-out in
  `src/workflows/run_benchmark/main.nf` (~L383). Present on `build/main` (config + `target/`).
- Sweep scripts: `scripts/run_benchmark/param_sweep/singler_params.yaml` +
  `run_test_singler_nebius.sh`.

## Optimization / tuning

**Grounded in the singler-py source** (`_train_single.py`, `_classify_single.py`,
fetched from `SingleR-inc/singler-py`), whose real defaults are:

`train_single(..., marker_method="classic", num_de=None, aggregate=False, ...)` and
`classify_single(..., quantile=0.8, use_fine_tune=True, fine_tune_threshold=0.05, ...)`.

### Tier 0 — the input, not a parameter
The **reference** dominates SingleR accuracy: which scRNA-seq atlas, how well its panel
overlaps the iST gene set, and **which annotation granularity** it is labelled at. The
mouse-brain test reference carries several nested label columns — `cell_type`,
`cell_type_level2`, `cell_type_level3`, `cell_type_level4` (coarse→fine). Selecting among
them is exactly what `--celltype_key` is *meant* to do (see Tier 1).

### Tier 1 — highest impact on output quality
- **`celltype_key` (annotation granularity)** — ALREADY EXPOSED but **UNWIRED** (see the
  "dead args" note). Coarser labels → higher per-class accuracy but less resolution; finer
  labels (`level3/4`) → more classes, harder correlation separation on a small iST panel.
  This is the one meaningful axis for a reference labeller and the basis of the prepared
  sweep — **but it only varies output after this one-line fix + a container rebuild**:

  ```python
  # L42, replace:
  ref_labels = sce_ref.get_column_data().column("cell_type")
  # with:
  ref_labels = sce_ref.get_column_data().column(par["celltype_key"])
  ```

- **`marker_method`** `{"classic","auc","cohens_d"}`, default `"classic"` — chooses how
  per-label markers are ranked. On small, low-count iST panels `cohens_d`/`auc` can pick
  more robust discriminative genes than classic pairwise log-FC. **Not exposed.**
- **`use_fine_tune`** (default `True`) / **`fine_tune_threshold`** (default `0.05`) —
  fine-tuning re-scores the top candidate labels using only their mutual markers; the
  single biggest accuracy lever between similar cell types. Currently left at the
  (good) default `True`. **Not exposed** (would need exposing to sweep the threshold or to
  measure the cost of turning it off).

### Tier 2 — quality/speed trade-offs
- **`num_de`** (default: auto for classic, else 10) — markers per pairwise comparison. **Not exposed.**
- **`quantile`** (default `0.8`) — quantile of the correlation distribution used for each
  label's score; lower is more robust to a few high-correlation outlier genes. **Not exposed.**
- **`aggregate`** (default `False`) — pseudo-bulk the reference for speed on large atlases,
  small accuracy change. **Not exposed.**

### Tier 3 — not exposed by the component (future work)
`marker_method`, `num_de`, `quantile`, `use_fine_tune`, `fine_tune_threshold`, `aggregate`
are all reachable in singler-py but **not surfaced** in `config.vsh.yaml`, and the script
passes none of them. Exposing any of these means: add the arg to `config.vsh.yaml`
`arguments:`, thread it into the `train_single`/`classify_single` call in `script.py`, then
`viash ns build` + rebuild the container (see the `check-component` skill) — a
**stale `build/main` container silently ignores a newly-added arg**, so a sweep over any of
these is NOT submittable until the rebuild lands.

### The prepared sweep (and its caveat)
`singler_params.yaml` sweeps **`celltype_key`** over `[cell_type_level2, cell_type_level3,
cell_type_level4]` (default `cell_type` covered by the default variant) → 4 variants total.
`celltype_key` already exists in `build/main`'s config, so the sweep **launches without a
rebuild**. HOWEVER, because of the dead-arg bug above, the deployed `build/main` script
hardcodes the `"cell_type"` column and ignores the value — so **until the one-line wiring
fix is applied and the container is rebuilt, the four variants produce identical output**
(the sweep is a no-op that measures nothing). Apply the fix + rebuild before using this
sweep to draw conclusions.

### Non-default audit
No exposed argument sits at a non-default (non-performance) value that deviates from an
upstream tool default: `celltype_key`/`labels_key` are data-key strings (no upstream tool
default to deviate from), and every singler-py algorithm knob is left at its library
default. Nothing was promoted onto the sweep axis by the non-default audit; `celltype_key`
is on the axis by design choice (Tier-1 granularity), not because it was mis-defaulted.
