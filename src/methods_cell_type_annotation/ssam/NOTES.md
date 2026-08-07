# ssam — cell-type annotation (SSAM)

## What this component is

Cell-type-annotation stage method (API `src/api/comp_method_cell_type_annotation.yaml`,
subtype `method_cell_type_annotation`). It labels each spatial cell with a cell type by the
**SSAM** approach: signature-based, cell-segmentation-free inference of cell types from an
mRNA-density map, then a per-cell majority vote over the transcripts assigned to that cell.

- Docs: https://ssam.readthedocs.io — Repo: https://github.com/HiDiHlabs/ssam
- Paper (config `references.doi`, verified correct): Park, Choi, Tiesmeyer et al.,
  "Cell segmentation-free inference of cell types from in situ transcriptomics data",
  *Nat Commun* 12, 3545 (2021), DOI `10.1038/s41467-021-23807-4` (PMID 34112806).
  No citation caveat — the config DOI is the SSAM method paper.

**Important implementation caveat.** Despite the name/DOI, this component does **not** run
the original `ssam` package. It calls `txsim.preprocessing.run_ssam` (txsim from
`theislab/txsim@dev`), which internally uses the **`plankton`/`planktonspace`**
re-implementation: `from plankton.utils import ssam`. The Docker setup installs
`planktonspace` + `matplotlib<3.9` on top of the spatialdata/txsim base. So the algorithm
is plankton's `ssam()`, with defaults that differ from the paper (see below).

## script.py — step by step

1. Assert `input_transcript_assignments` and `input_scrnaseq_reference` are provided
   (both required for this method).
2. Read `input_spatial_normalized_counts` (h5ad), the `transcripts` element of the
   `transcript_assignments.zarr`, and the SC reference; set `adata_sc.X` to its
   `layers["normalized"]`.
3. Subset the spatial AnnData to genes shared with the SC reference.
4. Call `tx.preprocessing.run_ssam(adata_sp, transcripts.compute(), adata_sc,
   um_p_px=par['um_per_pixel'], cell_id_col='cell_id', gene_col='feature_name',
   sc_ct_key=par['celltype_key'])`.
5. Copy `obs['ct_ssam']` -> `obs['cell_type']` (string) and write output.

**What `run_ssam` actually does** (txsim `preprocessing/_ctannotation.py::run_ssam`):
builds a `plankton.SpatialData` from the transcript gene column and `x*um_p_px`,
`y*um_p_px`; derives mean expression signatures per SC cell type; calls
`ssam(sdata, signatures=..., kernel_bandwidth=4, patch_length=1500, threshold_cor=0.2,
threshold_exp=0.1)` to produce a cell-type map; samples the map at every molecule; then
majority-votes a cell type per `cell_id`. `ct_ssam_cert` records the fraction of a cell's
spots that agree with the winning label.

Two live `# TODO`s in the script flag a known correctness risk: transcripts are passed in
**physical (µm) space**, not pixel space, so `um_p_px` scaling interacts with the fixed
kernel in a way the author suspects yields poor results ("ssam most likely outputs bad
results since the transcripts are provided in physical space instead of pixel space"). This
is exactly why `um_per_pixel` is the meaningful thing to sweep.

## Arguments

| Arg | Type | Config default | txsim/tool default | Maps to |
|-----|------|----------------|--------------------|---------|
| `--um_per_pixel` | double | **0.5** | `um_p_px=0.325` | scales transcript x/y before the density map (only exposed tuning knob) |
| `--celltype_key` | string | `cell_type` (from API) | `sc_ct_key='celltype'` | which SC `obs` column holds cell types (schema arg, not a tuning knob) |
| `--input_*` / `--output` | file | — | — | standard stage I/O |

Note: the config `um_per_pixel` default (0.5) **deviates** from txsim's `um_p_px` default
(0.325). The `# TODO` on the arg ("Should be able to infer this from transcripts") indicates
it is a placeholder rather than a tuned value.

## Setup / Docker

Merges `setup_spatialdata_partial.yaml` + `setup_txsim_partial.yaml` on
`openproblems/base_python:1`, then `pypi: [planktonspace, "matplotlib<3.9"]`. The
`matplotlib<3.9` pin is load-bearing (see commit `80c18d704 "ssam matplotlib"`): plankton's
plotting import breaks against matplotlib >= 3.9. No further Docker drama recorded.

## Wiring

- Registered as `celltype_annotation_methods` in the run_benchmark workflow config; the
  stage default is `tacco`. `main.nf` fans out one variant per enabled annotation method.
- Param sweep: `scripts/run_benchmark/param_sweep/ssam_params.yaml` +
  `run_test_ssam_nebius.sh` (CPU compute env, no `gpu` label; enables `tacco` + `ssam`).

## Optimization / tuning

**Tier 0 — input / coordinate space (biggest lever, not a sweepable arg).** Transcripts
arrive in µm (physical) space, contradicting SSAM's pixel-grid assumption. Fixing the
coordinate convention (or deriving the true µm/px from the transcript table, per the arg's
TODO) would likely matter more than any single-knob sweep.

**Tier 1 — the only exposed knob: `um_per_pixel`.** With the SSAM Gaussian
`kernel_bandwidth` fixed at 4 (µm-equivalent) inside txsim, `um_per_pixel` is the *only*
reachable control over the effective KDE smoothing scale: it rescales the point cloud
relative to the fixed kernel. Smaller -> coarser smoothing; larger -> finer. Swept over
`[0.1, 0.2, 0.325, 1.0]` (straddling the txsim default 0.325 and the shipped 0.5, which the
default variant already covers).

**Tier 3 — high-value knobs NOT exposed today (future work; NOT in this sweep).** These are
the true SSAM quality dials, but they are **hardcoded inside txsim's `run_ssam`** call to
`ssam(...)`, so the component cannot forward them. Exposing them means editing
`txsim/preprocessing/_ctannotation.py` (add params to the `run_ssam` signature + the `ssam()`
call), then adding matching `config.vsh.yaml` args and rebuilding the image
(`viash ns build` + container rebuild on `build/main` — see the `check-component` skill).
Until then they must NOT go in `ssam_params.yaml`:

- `kernel_bandwidth` (txsim hardcodes **4**; SSAM's own default is 2.5 µm) — the KDE
  bandwidth; the single biggest quality dial. Lower = sharper/noisier, higher = smoother.
- `threshold_cor` (hardcoded **0.2**) — minimum correlation between a pixel's local
  expression vector and a signature for a cell type to be called; SSAM typically uses ~0.6.
  Directly trades recall vs precision of assignments.
- `threshold_exp` (hardcoded **0.1**) — minimum total expression for a pixel to be
  classified at all (a foreground/vector-field threshold).
- `patch_length` (hardcoded **1500**) — tiling size; primarily a memory/speed knob (would
  stay fixed even if exposed).
- `no_ct_assigned_value` (`run_ssam` default `'None_sp'`) — label for unassigned cells;
  behavioural, not a quality dial.

## Risk points / gotchas

**`ZeroDivisionError` in `run_ssam` on datasets with < 20 cells (confirmed, upstream txsim
bug).** txsim `preprocessing/_ctannotation.py::run_ssam` throttles a progress print with

```python
n_cell_ids = len(adata_st.obs['cell_id'])
incr = n_cell_ids // 20          # == 0 when fewer than 20 cells reach annotation
for cell_id in adata_st.obs['cell_id']:
    if (count % incr) == 0:      # ZeroDivisionError: integer modulo by zero
```

When the spatial object entering annotation has fewer than 20 cells, `incr` floor-divides to
0 and `count % 0` raises on the first loop iteration. `incr` gates **only** a progress
`print` — it does not touch the annotation, so this is a pure robustness bug, independent of
`um_per_pixel` (every ssam variant fails identically on a too-small dataset). Fix is a
one-liner **in txsim** (not this repo): `incr = max(1, n_cell_ids // 20)`. Still unfixed on
`theislab/txsim@dev` HEAD as of 2026-08-06 (container ran `dev@0f62644b`, v0.1.2). Deploying
the fix requires patching txsim's `dev` and rebuilding the ssam container (`build_main`);
no change to this component's `config.vsh.yaml`/`script.py` is needed. See the iteration log
in the memory file for the sweep that surfaced it.
