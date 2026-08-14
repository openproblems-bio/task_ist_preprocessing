# baysor — transcript assignment (NOTES)

Baysor (Petukhov et al. 2021, Nat Biotechnol, DOI 10.1038/s41587-021-01044-w) re-segments
imaging-based spatial transcriptomics data by jointly modelling the likelihood of the local
transcriptional composition together with cell size/shape. In this component it is used as a
*transcript-assignment* method: a prior segmentation (the upstream `segmentation.zarr` label
map) is fed in as a per-transcript `cell_id`, and Baysor refines/assigns transcripts around it.

## How the component runs

- Wrapper: `script.py` runs Baysor **v0.7.1** through **sopa** (`sopa.segmentation.baysor`),
  NOT the Baysor CLI directly. CPU-only (no GPU).
- Prior: the script first assigns each transcript a `cell_id` by rasterizing the upstream
  segmentation label image at the transcript's (x, y), then hands that to sopa as
  `prior_shapes_key="cell_id"`. Baysor uses it as the prior segmentation.
- Patches: `sopa.make_transcript_patches(patch_width=2000, patch_overlap=50)` is **hardcoded** —
  a memory/parallelism knob, not a quality knob. Leave fixed.
- Threads: `JULIA_NUM_THREADS = max(cpus-2, 1)` — resource knob, leave fixed.
- The six scientific knobs are all forwarded from `par[...]` into sopa's Baysor config dict:
  - `data.force_2d`, `data.min_molecules_per_cell`
  - `segmentation.scale`, `segmentation.scale_std`, `segmentation.n_clusters`,
    `segmentation.prior_segmentation_confidence`
  All six are already exposed as `config.vsh.yaml` `arguments:` — nothing is hardcoded that
  needs exposing.

## Coordinate units

Global coordinate space here is ~1 µm/px, so `scale` (Baysor's "typical cell radius, in the
same units as x/y") is effectively **µm**. `scale = -1.0` is the sentinel for "estimate the
radius from `min_molecules_per_cell`" (Baysor: if `scale <= 0` it is estimated). sopa's own
example config instead sets an explicit `scale: 6.25` µm.

---

## Optimization / tuning

Goal: assignment quality (correct, well-shaped cells; good transcript-to-cell agreement with
the reference), CPU time secondary. Config defaults here are NOT all Baysor defaults — two were
deliberately moved (see audit), so tuning mostly means walking those back toward the tool
default while also probing the two dominant size knobs (`scale`, `min_molecules_per_cell`).

**Tool-default audit (config default → Baysor/sopa true default):**

| arg | config default | Baysor CLI default | sopa example | verdict |
|-----|---------------|--------------------|--------------|---------|
| `force_2d` | `true` | (2D mode) | `true` | matches sopa; 2D benchmark data → **fixed** (boolean flip to 3D is meaningless here) |
| `min_molecules_per_cell` | `50` | none ("important param") | `10` | **DEVIATES** → swept |
| `scale` | `-1.0` (auto) | none (required or est.) | `6.25` µm | deliberate auto mode → swept (explicit radii) |
| `scale_std` | `"25%"` | `"25%"` | `"25%"` | matches → Tier-2 axis only |
| `n_clusters` | `4` | `4` | (unset) | matches → Tier-2 axis only |
| `prior_segmentation_confidence` | `0.8` | `0.2` | `0.5` | **DEVIATES** → swept |

Structural args `transcripts_key`, `coordinate_system` are not tuning knobs (not swept).

### Tier 0 — input (not a parameter)
The prior segmentation quality itself (upstream `segmentation.zarr`) and coordinate scaling are
the biggest levers. Fixed by the pipeline; called out for context.

### Tier 1 — highest impact (all in the sweep)
- **`scale`** — target cell radius (µm here). The single most important Baysor parameter.
  Default `-1.0` = estimate from `min_molecules_per_cell`. Sweep explicit radii
  **[5.0, 7.5, 10.0]** µm, grounded in typical mammalian cell radii (~4–10 µm) and straddling
  sopa's `6.25`. Tests explicit-radius vs auto-estimate.
- **`prior_segmentation_confidence`** — trust in the prior segmentation vs re-segmenting from
  scratch (∈[0,1]). **Deviates**: config `0.8` (Baysor: ">0.7 only when the prior is highly
  reliable") vs Baysor default `0.2`. Sweep **[0.2, 0.5]** — Baysor default and the sopa-example
  value, walking back from the trust-the-prior extreme toward the tool default.
- **`min_molecules_per_cell`** — minimum molecules to call a cell; also drives the auto scale
  estimate. **Deviates**: config `50` vs sopa example `10`. Sweep **[10, 25]** (toward the
  smaller sopa value; more permissive cell calling + smaller estimated radius).

### Tier 2 — companion quality knobs (at Baysor default; probed lightly)
- **`scale_std`** — prior std of cell radius, relative to `scale`. Default `"25%"`. Sweep
  **["50%"]** (wider size tolerance).
- **`n_clusters`** — number of molecule clusters (≈ major cell types) used in the neighbourhood
  composition prior. Default `4`; Baysor recommends 3–15. Sweep **[6, 8]**.

### Tier 3 — not exposed (none needed)
No high-value Baysor knob is hidden in the script: all six scientific parameters are already
config args. Other Baysor knobs (`min_molecules_per_segment`, `confidence_nn_id`,
`estimate_scale_from_centers`, `nuclei_genes`/`cyto_genes`) are lower-value and left unexposed.

### Status
**SUBMITTABLE AS-IS** — every swept arg already ships in `config.vsh.yaml`, so the
`--revision build/main` container accepts them with **no rebuild** and no config/script edits.

Sweep = 1 default + 3 (`scale`) + 2 (`prior_segmentation_confidence`) + 2
(`min_molecules_per_cell`) + 1 (`scale_std`) + 2 (`n_clusters`) = **11 variants**.

Fixed resource/perf knobs (never swept): `make_transcript_patches` `patch_width=2000` /
`patch_overlap=50`, `JULIA_NUM_THREADS`, and `force_2d` (2D data).
