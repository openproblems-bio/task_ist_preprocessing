# Tangram cell-type annotation — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how it
works, which of Tangram's many knobs are actually wired up, why the two shipped
defaults were chosen, and where to tune for quality vs. time.

Links: docs https://tangram-sc.readthedocs.io · repo
https://github.com/broadinstitute/Tangram.

**Citation:** the config's `references.doi` is `10.1038/s41592-021-01264-7` —
Biancalani et al., "Deep learning and alignment of spatially resolved single-cell
transcriptomes with Tangram", *Nature Methods* 18, 1352–1362 (2021), PMID
34711971. Verified via PubMed: this **is** the method-specific paper (not a generic
framework reference), so no citation caveat is needed.

## What this component is

One of the methods at the **cell-type-annotation** stage of the iST preprocessing
benchmark (`src/methods_cell_type_annotation/`). It implements the standard API
`src/api/comp_method_cell_type_annotation.yaml`:

- in: `input_spatial_normalized_counts` (`.h5ad`, cell × gene with a `normalized`
  layer) + `input_scrnaseq_reference` (`.h5ad`, reference atlas with a `normalized`
  layer and a `celltype_key` obs column, default `cell_type`).
- out: the spatial AnnData with a per-cell label written to `obs[celltype_key]`.

It is a **GPU deep-learning method**: Tangram (torch) learns a soft mapping matrix
between the reference and the spatial cells by gradient descent, then projects the
reference labels through that mapping. Nextflow directives are
`[ midtime, midcpu, midmem, gpuh100 ]` — the `gpuh100` label is what schedules it
onto a GPU node.

## script.py — step by step

1. **Device pick** (`:19-23`) — `"cuda:0"` if `torch.cuda.is_available()` else
   `"cpu"`. Tangram runs on CPU but is much slower there.
2. **Read inputs** (`:30-31`) — spatial + reference AnnData.
3. **Use the log1p-normalized layer as `X`** (`:33-35`) — both adatas' `X` is
   replaced by their `normalized` layer before mapping. (That layer is produced
   upstream by the normalization stage / the SC process workflow, not here.)
   `adata_sp_orig` is stashed so the final output carries the original layers, not
   Tangram's intermediate additions.
4. **Markers = the full spatial panel** (`:39-40`) — every `adata_sp.var_name` is
   used as a training gene. iST panels are small (hundreds of genes) so "all genes"
   is reasonable; there is no marker-selection step.
5. **`tg.pp_adatas`** (`:44`) — intersects genes across the two adatas, drops
   all-zero genes, and stores density priors. Mutates both adatas in place.
6. **Build `map_kwargs` and map** (`:56-66`) — calls `tg.map_cells_to_space` with
   `mode=par['mode']`, `num_epochs=par['num_epochs']`, `density_prior='uniform'`
   (hardcoded), `device`. In `clusters` mode it also passes
   `cluster_label=par['celltype_key']` (`:64-65`) — required so Tangram averages
   reference cells within each cluster.
7. **Project labels** (`:69-73`) — `tg.project_cell_annotations` writes a
   spot × cell-type score frame to `adata_sp.obsm['tangram_ct_pred']`.
8. **Argmax → label** (`:76-80`) — restore `adata_sp_orig`, then set
   `obs[celltype_key] = tangram_ct_pred.idxmax(axis=1)` (hard label = highest-scoring
   type). The commented block (`:83-87`) would additionally store a per-cell
   confidence score and the full score matrix — currently disabled.
9. **Write** (`:90`).

## Arguments

Only **two** knobs are exposed; everything else in `map_cells_to_space` is left at
the Tangram default or hardcoded in the script.

| config arg | default | Tangram default | forwarded to | note |
|---|---|---|---|---|
| `--mode` | `clusters` | `cells` | `map_cells_to_space(mode=…)` | **deviates** from tool default; see below |
| `--num_epochs` | `1000` | `1000` | `map_cells_to_space(num_epochs=…)` | matches tool default |

**Why `mode: clusters` (a deliberate non-default).** In `cells` mode Tangram learns
an `(n_sc_cells × n_spatial_cells)` mapping matrix on the GPU, which exhausts VRAM
for large references (the config comment cites the 101k-cell MPII lung reference).
`clusters` mode maps per-cell-type *clusters* instead, shrinking the matrix by
~cells-per-type so it fits comfortably. It is also the natural mode for **label
projection** (we only want the cell-type label, not an individual-cell
correspondence). So the deviation is a VRAM/scale choice — but note it *does* change
the output (cluster-level vs. cell-level mapping), so it is treated as a tunable
axis, not a silent baseline (see "Optimization / tuning").

**Hardcoded (not a config arg):** `density_prior='uniform'` (`:62`). Tangram's tool
default is `'rna_count_based'`. The docstring guidance (`:48-50`) is: use `'uniform'`
when the spatial voxels are at **single-cell resolution** (MERFISH/Xenium-style),
and `'rna_count_based'` when density is expected to track RNA counts. Since iST here
is single-cell-resolution, `'uniform'` is the appropriate choice — but it is a
deviation from the tool default and is **not exposed**, so it can only be swept
after being added to the config (see Tier 3).

## config.vsh.yaml — the Docker setup

- Base image `openproblems/base_python:1` + `pypi: [tangram-sc]`.
- The config carries a **history comment**: `openproblems/base_pytorch_nvidia:1`
  (the CUDA/torch base used by other GPU methods) was tried first and "leads to
  dependency issues", so it was reverted to the plain python base and lets
  `tangram-sc` pull its own torch. A `TODO` notes a different pytorch+CUDA base
  might be worth trying. If you touch the image, that dependency clash is the thing
  to re-verify.
- `native` engine is also declared (for local/CI without Docker).

## Optimization / tuning

Two framing facts before tuning:
1. Unlike the segmentation methods, the tangram defaults are **not** aggressively
   speed-floored — `num_epochs=1000` is Tangram's own recommended default. The one
   real deviation is `mode` (chosen for VRAM/scale, see above).
2. All ranges below are grounded in the Tangram docs signature (verified against
   `tangram/mapping_utils.py`: `map_cells_to_space(mode="cells",
   learning_rate=0.1, num_epochs=1000, lambda_g1=1, lambda_d=0,
   density_prior='rna_count_based', …)`).

**Tier 0 — the input, not a parameter.** The biggest lever is the *reference*: which
SC atlas, how well its cell types match the tissue, and the shared gene panel.
`pp_adatas` restricts training to genes present in both — a poorly matched reference
caps achievable accuracy regardless of any knob. Also note markers = the entire
spatial panel (no marker curation).

**Tier 1 — highest impact on quality (quality vs. time): `num_epochs`**
*(exposed; default 1000 = tool default).* Tangram optimizes the mapping by gradient
descent for `num_epochs` steps; more epochs = better-converged mapping (to a point),
fewer = faster but under-converged and noisier labels. This is the primary
accuracy↔runtime dial. Sweep below and around 1000 (e.g. 100 / 500 / 2000 / 3000) to
find the knee — where extra epochs stop improving the argmax labels.

**Tier 2 — `mode`** *(exposed; default `clusters`, tool default `cells`).* Not just a
resource knob — it changes the mapping granularity and therefore the labels.
`cells` maps individual reference cells (finer, can be more accurate on **small**
references) but its matrix scales with `n_sc_cells` and OOMs on large ones; the
config picked `clusters` for that reason. The one meaningful non-default value to
sweep is `cells` (feasible on the small test reference; expect OOM on large
references such as MPII lung on limited VRAM). Tangram also offers `constrained`,
but that mode is designed for **cell filtering** and needs `target_count` +
`lambda_count`/`lambda_f_reg` to be meaningful, so it is not a drop-in axis for pure
label projection — left out of the sweep.

**Tier 3 — high-value knobs NOT exposed today (future work; would need config
expose + container rebuild, so NOT in the submittable sweep):**
- **`density_prior`** — hardcoded `'uniform'`; tool default `'rna_count_based'`.
  Directly shapes how mass is distributed across spatial cells. Expose it (choices
  `uniform` / `rna_count_based`) to test whether the single-cell-resolution
  assumption actually helps on each dataset. Highest-value un-exposed knob.
- **`learning_rate`** *(tool default 0.1)* — step size of the Adam optimizer;
  interacts with `num_epochs` (lower LR needs more epochs). Worth a small sweep
  (0.05 / 0.1 / 0.5) once epochs are set.
- **`lambda_g1` / `lambda_d` / `lambda_r`** *(1 / 0 / 0)* — the loss weights
  (gene-expression term, density term, entropy regularizer). Advanced; only worth
  touching after epochs + density_prior are settled.

To sweep any Tier-3 knob you must add it to `config.vsh.yaml`'s `arguments:` (type +
default = the Tangram default), thread it into the `map_kwargs` dict in `script.py`,
then `viash ns build` + rebuild the container (see the `check-component` skill). A
sweep launched against a stale `build/main` container would silently ignore a
freshly-added arg — which is exactly why these stay out of the current
submittable-now sweep.

## Risk points / gotchas

- **`cells` mode OOMs on large references.** The whole reason `clusters` is the
  default. Sweeping `mode: cells` is safe on the small test reference but expect
  VRAM blow-ups on 100k-cell references at limited VRAM.
- **`density_prior` is hardcoded, not a config arg.** Changing it requires editing
  `script.py`, not the sweep file.
- **Only a hard argmax label is emitted** (`:80`); the per-cell confidence score
  block is commented out (`:83-87`). Downstream metrics see labels only, no scores.
- **GPU base-image fragility.** The plain `base_python:1` + `tangram-sc` combo is
  deliberate; `base_pytorch_nvidia:1` caused dependency issues (see config comment).
- **No end-to-end tuning run has been done yet** — this NOTES.md documents the setup
  and the sweep design; results are not yet in.

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: dependency
  `methods_cell_type_annotation/tangram` (`:181`); part of the annotation default
  string `ssam:tacco:moscot:mapmycells:tangram:singler:rctd` (`:101`).
- `src/workflows/run_benchmark/main.nf:382`: imported in the annotation fan-out.
- Run scripts: enabled (uncommented) in `run_gpu_nebius.sh`; commented in the plain
  `run_test_nebius.sh` / `run_test_local.sh` (GPU + time cost).
- **Parameter sweep:** `scripts/run_benchmark/param_sweep/tangram_params.yaml`
  (committed default + sweep) and `run_test_tangram_nebius.sh` (Nebius GPU launch,
  reads the params file from GitHub raw). The stage default `tacco` is kept alongside
  `tangram` so the run has a baseline to compare against.

## References

- Docs: https://tangram-sc.readthedocs.io (the
  `tangram.mapping_utils.map_cells_to_space` page has the full parameter reference).
- Repo: https://github.com/broadinstitute/Tangram
- **Paper (the DOI in `config.vsh.yaml`, verified correct):** Biancalani T. et al.,
  "Deep learning and alignment of spatially resolved single-cell transcriptomes with
  Tangram", *Nature Methods* 18, 1352–1362 (2021), DOI 10.1038/s41592-021-01264-7,
  PMID 34711971.
