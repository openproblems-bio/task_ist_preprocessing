# TACCO (tacco) — cell-type-annotation method notes

## What it is
TACCO ("Transfer of Annotations to Cells and their COmbinations") transfers cell-type
labels from an annotated scRNAseq **reference** onto the spatial cells by modelling the
data as **continuous mixtures** and solving an **optimal-transport (OT)** problem, wrapped
with optional boosters/normalizers. It matches or beats specialised annotators while being
up to ~10× faster.

- Paper: Mages, Moriel, Avraham-Davidi, Murray, Watter, Chen, Rozenblatt-Rosen,
  Klughammer, Regev, Nitzan. **"TACCO unifies annotation transfer and decomposition of
  cell identities for single-cell and spatial omics."** *Nat Biotechnol* 41, 1465–1473
  (2023). DOI [10.1038/s41587-023-01657-3](https://doi.org/10.1038/s41587-023-01657-3).
- Docs: <https://simonwm.github.io/tacco/> · Repo: <https://github.com/simonwm/tacco>
- **Citation is correct** — the config `references.doi` already points at the method's own
  paper (not a generic framework paper), so no citation caveat.

## How the component works (`script.py`)
1. Reads `input_spatial_normalized_counts` (spatial) + `input_scrnaseq_reference` (SC).
2. **Switches both to raw counts** (`adata.X = adata.layers['counts']`) — TACCO wants raw
   counts, not the normalized layer. This is a Tier-0 input fact: despite the arg name
   `*_normalized_counts`, tacco annotates on the **counts** layer.
3. Calls `tacco.tl.annotate(adata=adata_sp, reference=adata_sc, annotation_key=celltype_key, **knobs)`.
   `annotate` (with `result_key=None`) returns an `n_obs × n_celltypes` **proportion**
   DataFrame.
4. Takes the **argmax** per cell → hard label written to `adata_sp.obs[celltype_key]`.
5. Writes `output`.

`celltype_key` (default `cell_type`) is the reference `.obs` column holding the labels to
transfer, and also the output column name.

## Optimization / tuning
Before this pass the component forwarded **only** `annotation_key`; every algorithmic knob
ran at the tacco library default. All knobs below were **newly exposed** here
(`config.vsh.yaml` `arguments:` + threaded into the `annotate(**knobs)` call in
`script.py`). Each defaults to tacco's own default, so an unset value == stock behaviour.

Full `tacco.tl.annotate` default signature (from the docs): `method='OT'`,
`multi_center=None`, `platform_iterations=None`, `bisections=None`, `bisection_divisor=3`,
`multi_center_amplitudes=True`, plus gene/QC filters (`min_log2foldchange`, `n_hvg`,
`remove_mito`, `max_annotation`, …) left at defaults.

- **Tier 0 — input (not a knob).** Raw **counts** are fed (see above), and the reference is
  the panel-subset SC. Biggest lever is really reference quality/granularity (`celltype_key`
  level), but that is a data choice, not swept here (kept at `cell_type`, matching moscot/
  tangram; only singler sweeps the label level).
- **Tier 1 — highest impact on annotation quality:**
  - **`multi_center`** (default None = off). k-means-splits each reference type into
    sub-centers so the OT plan matches within-type substates instead of one mean profile —
    the paper's headline accuracy mechanism. Tutorials use ~10. **Swept `[2, 5, 10, 20]`.**
  - **`platform_iterations`** (default None = off). Corrects systematic gene-capture
    differences between the scRNAseq reference and the sparse iST panel before annotating —
    directly relevant to imaging ST. `0` = one normalization pass, `>0` iterates.
    **Swept `[0, 2, 3]`.**
  - **`method`** (default `'OT'`). The annotation wrapper. OT is the paper's core model;
    `nnls` and `projection` are dependency-free baselines present in this image.
    **Swept `[nnls, projection]`.** RCTD/SingleR/tangram/novosparc/WOT backends need extra
    packages the container does **not** install → NOT swept (would just fail).
- **Tier 1/2 — quality at compute cost:**
  - **`bisections`** (default None = off). Runs the annotator recursively, each round
    removing a reconstructed fraction of counts, to resolve mixed/low-count cells. Tutorials
    use `4` (with `bisection_divisor=3`). **Swept `[2, 3, 4]`.**
- **Tier 2/3 — exposed but NOT swept:**
  - **`bisection_divisor`** (default 3, exposed): only has an effect when `bisections > 0`;
    left at the tacco default and gated on `bisections` in `script.py`.
  - Gene-selection / QC filters (`min_log2foldchange`, `n_hvg`, `remove_mito`,
    `max_annotation`, `multi_center_amplitudes`) — left at tacco defaults; candidates for a
    deeper sweep later.

### Sweep files
- `scripts/run_benchmark/param_sweep/tacco_params.yaml` — star-around-default,
  **13 variants** (1 default + 4 + 3 + 3 + 2).
- `scripts/run_benchmark/param_sweep/run_test_tacco_nebius.sh` — CPU nebius launcher (no
  `gpu` label), reads the params yaml from its raw-GitHub URL on the current branch.

### ⚠ Deployment gate (read before launching)
Every swept knob is a **newly added config arg**, so the deployed `build/main` tacco
container does **not** have them yet. The sweep only runs after:
1. `viash ns build` (regenerate the target) — locally, to test; and
2. **merge to `main` → CI rebuilds the tacco container**; verify with the `check-component`
   skill (`check-component tacco` → fresh) before launching.
Until then, non-default variants error with `unknown option '--multi_center'` (only the
default variant, stock OT, would run). `run_sweep`'s build gate enforces this.

The params yaml must also be **committed and pushed** to the launch branch (the launcher
reads it from `raw.githubusercontent.com/.../<branch>/...`), same as the other sweeps.
