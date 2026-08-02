# moscot — cell-type annotation (MOSCOT / optimal transport)

## What this component is

Cell-type-annotation stage method (API `src/api/comp_method_cell_type_annotation.yaml`,
subtype `method_cell_type_annotation`). It labels each spatial cell by **optimal-transport
mapping** of an scRNA-seq reference onto the spatial cells and then transferring the
reference's cell-type labels through the learned transport plan.

- Docs: https://moscot.readthedocs.io — Repo: https://github.com/theislab/moscot
- Paper (config `references.doi`, **verified correct**): Klein, Palla, Lange, Klein et al.,
  "Mapping cells through time and space with moscot", *Nature* (2025),
  DOI `10.1038/s41586-024-08453-2`. This is the moscot method paper — **no citation caveat**.
- **GPU-only** in practice: the engine installs `jax[cuda12]` + `moscot` + `flax` + `diffrax`
  on `openproblems/base_pytorch_nvidia:1.1.0`, and the Nextflow directives carry `gpuh100`.

Under the hood it uses `moscot.problems.space.MappingProblem`, a **fused Gromov-Wasserstein
(FGW)** problem solved with OTT-jax: a quadratic (Gromov-Wasserstein) term matches
intra-domain structure (the SC feature space vs the spatial coordinates) and a linear term
matches shared-gene expression across the two domains; `alpha` interpolates between them.

## script.py — step by step

1. **`:13-15` — pop `LD_LIBRARY_PATH` before importing jax.** JAX's `jax[cuda12]` wheel ships
   its own CUDA/cuDNN; a system `LD_LIBRARY_PATH` pointing at an older cuDNN otherwise wins and
   crashes with "Loaded runtime CuDNN library: 9.1.0 but source was compiled with: 9.8.0".
   **Load-bearing** — do not remove.
2. `:18-25` — jax GPU sanity prints (version, backend, devices).
3. `:58-59` — read `input_scrnaseq_reference` and `input_spatial_normalized_counts` (h5ad).
4. **`:62-66` — SMALL-DATA GUARD (the key gotcha).** If `adata_sp.n_obs < 10000` it
   **overrides** `par['rank'] = -1` (full rank) and `par['tau'] = 1.0` (balanced OT),
   regardless of what was passed in. See "Risk points".
5. `:69-71` — assert a `"normalized"` layer exists in both AnnDatas and `centroid_x`/`centroid_y`
   exist in the spatial `obs`.
6. `:74-76` — set `X = layers["normalized"]` for both; build `adata_sp.obsm["spatial"]` from the
   centroid columns.
7. `:78` — `sc.pp.pca(adata_sc, n_comps=30)` — the SC quadratic cost is built from a **hardcoded
   30-dim PCA**, not raw genes.
8. `:81-86` — `MappingProblem(adata_sc, adata_sp).prepare(sc_attr={obsm, X_pca},
   xy_callback="local-pca")`: SC intra-domain cost = the 30-dim PCA; spatial intra-domain cost =
   a local-PCA callback over the spatial coordinates; the linear (fused) term links shared genes.
9. `:92-99` — `mp.solve(alpha=, epsilon=, tau_a=tau, tau_b=tau, rank=, batch_size=)`.
10. `:102-108` — `mp.annotation_mapping(mapping_mode=, annotation_label=celltype_key,
    source="src", forward=False)`; the returned per-cell label is written to
    `adata_sp.obs[celltype_key]`. (`forward=False`/`source="src"` are hardcoded.)
11. `:111` — write the annotated spatial AnnData.

## Arguments

| Arg | Type | Config default | moscot tool default | Effective on the <10k test data | Maps to |
|-----|------|----------------|---------------------|---------------------------------|---------|
| `--alpha` | double | **0.8** | **0.5** | used as-is | `solve(alpha=)` — FGW quadratic↔linear weight |
| `--epsilon` | double | 0.01 | 0.01 | used as-is | `solve(epsilon=)` — entropic regularization |
| `--tau` | double | **0.3** | **1.0** | **forced to 1.0** (`:65-66`) | `solve(tau_a=tau_b=)` — marginal relaxation / unbalancedness |
| `--rank` | integer | **500** | **-1** | **forced to -1** (`:63-64`) | `solve(rank=)` — low-rank OT (−1 = full rank) |
| `--batch_size` | integer | 1024 | `None` | used (but n_obs<batch ⇒ no effect) | `solve(batch_size=)` — online GW-cost batching (memory knob) |
| `--mapping_mode` | string (`sum`/`max`) | `max` | *(required arg, no tool default)* | used as-is | `annotation_mapping(mapping_mode=)` — label read-out from the plan |
| `--celltype_key` | string | `cell_type` (from API) | — | — | which reference `obs` column supplies labels (schema arg) |

**Defaults history / why they deviate.** Originally `tau=1.0, rank=-1`. `rank` was raised to a
5000 default (`d3f45813f`, "for 300k cells") with the `<10k → rank=-1` guard added at the same
time; then `Adjust moscot parameters for memory efficiency` (`b10ce0563`) lowered the defaults to
`tau=0.3, rank=500` and extended the guard to also clamp `tau=1.0` on small data. So the config's
`tau=0.3`/`rank=500` are **large-dataset memory settings** that never take effect on small inputs.

## Setup / Docker

Merges nothing custom — the engine is `openproblems/base_pytorch_nvidia:1.1.0` plus
`pypi: [numpy, "jax[cuda12]", anndata, scanpy, moscot, flax, diffrax]`. History worth knowing:
`flax` (`997a9114b`) and `diffrax` (`d11af7b58`) are **explicit deps moscot needs at import**;
GPU container was switched in `702ba3fd2` / `ecb67b533`. The `LD_LIBRARY_PATH` unset in
`script.py` (not Docker) is the load-bearing cuDNN fix.

## Wiring

- Registered as a `celltype_annotation_methods` option in the run_benchmark workflow config; the
  stage default is `tacco`. `main.nf` fans out one variant per enabled annotation method.
- Param sweep: `scripts/run_benchmark/param_sweep/moscot_params.yaml` +
  `run_test_moscot_nebius.sh` (GPU compute env + `gpu` label; enables `tacco` + `moscot`).

## Optimization / tuning

Two framing facts first:
1. All tool defaults below were read from moscot's own source/docs
   (`MappingProblem.solve`: `alpha=0.5, epsilon=1e-2, tau_a=tau_b=1.0, rank=-1,
   scale_cost="mean", batch_size=None`).
2. **The test resource is tiny (306 spatial cells), so the `<10k` guard is always active on it.**
   That directly shapes which knobs are worth sweeping *on the test run*.

**Tier 0 — the input, not a parameter (biggest lever).** The reference atlas (how well its cell
types match the tissue) and the shared gene panel cap achievable accuracy. Also hardcoded and
therefore Tier-0/Tier-3: the SC representation fed to the quadratic cost is a **30-dim PCA**
(`sc.pp.pca(n_comps=30)`) and the spatial cost is a `local-pca` callback — changing either changes
*what structure* OT matches, independent of any exposed knob.

**Tier 1 — highest impact on quality, and EFFECTIVE on the test data (swept):**
- **`alpha`** *(exposed; config 0.8, tool default 0.5).* FGW interpolation: `alpha→1` pure
  Gromov-Wasserstein (intra-domain structure), `alpha→0` pure linear (shared-gene expression).
  The config's 0.8 is structure-heavy and **deviates from the tool default** ⇒ promoted to a
  sweep axis. Range `[0.5, 0.7, 0.9]` straddles the tool default up toward pure-GW; 0.8 omitted
  (the default variant covers it). (Author `# TODO` in config suggested `[0.7, 0.8, 0.9]`.)
- **`epsilon`** *(exposed; config 0.01 = tool default).* Entropic regularization / plan
  stochasticity: higher = blurrier map, faster & more stable convergence; lower = sharper,
  more confident mapping but harder to converge. Not a deviation, but the single most important
  OT knob, so swept as a legitimate quality axis: `[0.001, 0.005, 0.05, 0.1]` (an order of
  magnitude either side of 0.01).

**Tier 1 for LARGE datasets, but NEUTRALIZED on the test data (documented, NOT swept):**
- **`tau`** (`tau_a=tau_b`, unbalancedness/marginal relaxation) and **`rank`** (low-rank OT).
  Both deviate from tool defaults (`tau` 0.3 vs 1.0, `rank` 500 vs −1), which normally makes
  them mandatory sweep axes. **But `script.py:62-66` clamps `tau=1.0`/`rank=-1` on any input with
  `<10000` cells**, and the test resource has 306 — so on the submittable test run they *cannot
  change the output* (they behave like the skill's fixed "can't change output" exception).
  Sweeping them here would only add identical no-op variants. They are the highest-value knobs on
  a **real (>10k-cell) dataset**: for a real-data sweep, `tau ∈ [0.1, 0.2, 0.3, 1.0]` (the config
  `# TODO` notes it "seems only to work with tau=1 on our data") and `rank ∈ [500, 1000, 2000,
  5000, -1]` (the `# TODO` scales rank with cell count, ~5000 for 300k cells).

**Tier 2 — label read-out (swept):**
- **`mapping_mode`** *(exposed; config `max`, moscot has no tool default — required arg).*
  `max` = label of the single highest-mass source cell; `sum` = aggregate plan mass per cell type
  then argmax (uses the full coupling, smoother). One meaningful non-default value: `[sum]`.

**Tier 3 — high-value knobs NOT exposed today (future work; would need config expose + script
wiring + a container rebuild, so NOT in the submittable sweep):**
- **`n_comps`** — the SC PCA dimensionality is hardcoded at 30 (`:78`); it sets the resolution of
  the quadratic cost. Highest-value un-exposed knob.
- **`scale_cost`** *(tool default "mean")* — hardcoded (unset ⇒ tool default); how the cost
  matrices are normalized before the solve, affects the effective `epsilon`.
- **`initializer`** — for low-rank solves (`rank>0`) the initializer (e.g. `"rank2"`/`"k-means"`)
  materially changes convergence; only relevant once `rank` is swept on a large dataset.
- **`scale_by_marginals`** *(annotation_mapping, tool default True)* and **`forward`/`source`**
  (hardcoded `forward=False`, `source="src"`) — the label-projection direction/weighting.
- **`threshold` / `max_iterations`** *(convergence: tool `threshold=1e-3`)* — Sinkhorn/LR-GW
  convergence controls (quality↔runtime).

To sweep any Tier-3 knob: add it to `config.vsh.yaml` `arguments:` (default = the moscot default),
thread it into the `solve`/`prepare`/`annotation_mapping` call in `script.py`, then
`viash ns build` + rebuild the container (see `check-component`). A sweep launched against a stale
`build/main` container silently ignores a freshly-added arg — which is why these stay out of the
submittable-now sweep.

## Risk points / gotchas

- **The `<10k` small-data guard silently overrides `tau` and `rank`** (`:62-66`). On the 306-cell
  test resource this means those two args do nothing — any sweep over them is a no-op. This is
  *intended* behavior (unbalanced/low-rank OT is only for large references), not a bug.
- **GPU-only + not installable/verifiable locally.** moscot/jax[cuda12] is not in the local env;
  all tool defaults and signatures in this NOTES were read from moscot's upstream source and docs,
  not executed here. The sweep has **not** been run end-to-end yet.
- **A hard argmax label is emitted** — the transport plan / per-cell confidence is not exported;
  downstream metrics see labels only.
- **`batch_size` is a pure memory knob** (online GW-cost batching to bound GPU memory on large
  references); it does not change the output and is held fixed (`n_obs=306 < 1024` ⇒ no batching
  on the test data anyway).
