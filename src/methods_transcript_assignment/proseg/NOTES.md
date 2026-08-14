# proseg — NOTES

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the `proseg`
transcript-assignment method. Committed with the code.

## What this component is

- **Stage / API:** `methods_transcript_assignment` (merges
  `/src/api/comp_method_transcript_assignment.yaml`). Takes `raw_ist.zarr` + a prior
  `segmentation.zarr` (nuclear/cell label image) and reassigns every transcript to a cell,
  writing a transcripts-with-`cell_id` zarr.
- **Underlying tool:** **Proseg** ("Probabilistic Segmentation", Rust CLI, installed via
  `cargo install proseg`), driven through **sopa** (`sopa.segmentation.proseg`). proseg runs a
  probabilistic sampler inspired by Cellular-Potts models over a voxel lattice, using the prior
  segmentation as a nuclear prior, and jointly infers cell morphology + transcript assignment
  from the spatial distribution of transcripts. CPU (multi-threaded); no GPU.
- **How args reach the tool:** `script.py` builds a `command_suffix` string of `--flag value`
  pairs and passes it to `sopa.segmentation.proseg(..., command_line_suffix=...)`, which appends
  it to the proseg CLI. `--nthreads` is computed from `meta['cpus']` (`cpus-2`, min 1).
- **Links:** docs/repo https://github.com/dcjones/proseg .

## Citation

`references.doi: 10.1038/s41592-025-02697-0` — Jones et al., 2025, *Nature Methods*, the proseg
method paper. Method-specific, not a generic framework citation. No caveat needed.

## Optimization / tuning

**Verified upstream defaults** (proseg `src/main.rs`, `main` branch — the clap arg defaults and
the `DEFAULT_*` consts at L38-41): `ncomponents=10`, `voxel_layers=4`, `burnin_samples=200`,
`samples=200`, `cell_compactness=0.04`, `nuclear_reassignment_prob=0.2`,
`diffusion_probability=0.2`, `diffusion_sigma_near=1.0`, `diffusion_sigma_far=4.0`,
`voxel_size=1.0` (`DEFAULT_VOXEL_SIZE`), `burnin_voxel_size=2.0` (`DEFAULT_BURNIN_VOXEL_SIZE`).

**Newly exposed knobs (this branch).** Six proseg model knobs were HARDCODED in `script.py`'s
`command_suffix` and commented "these should not need to be changed". They are now proper
config `arguments:` (threaded into `command_suffix` as `{par['<arg>']}`, matching the existing
`--voxel-size {par['voxel_size']}` pattern), with defaults set to proseg's own CLI defaults:

| config arg (snake_case)      | proseg CLI flag              | new default | was hardcoded |
|------------------------------|------------------------------|-------------|---------------|
| `cell_compactness`           | `--cell-compactness`         | 0.04        | **0.03**      |
| `nuclear_reassignment_prob`  | `--nuclear-reassignment-prob`| 0.2         | 0.2           |
| `diffusion_probability`      | `--diffusion-probability`    | 0.2         | 0.2           |
| `ncomponents`                | `--ncomponents`              | 10          | 10            |
| `diffusion_sigma_far`        | `--diffusion-sigma-far`      | 4.0         | 4             |
| `diffusion_sigma_near`       | `--diffusion-sigma-near`     | 1.0         | 1             |

Because these are new config arguments, **the sweep NEEDS REBUILD**: a `build/main` Nebius run
uses a container that predates these args and will reject them until the component is
regenerated (`viash ns build`) and its container rebuilt (see `check-component`).

**Non-default audit result:**
- The five args exposed before this branch (`voxel_layers`, `samples`, `burnin_samples`,
  `voxel_size`, `burnin_voxel_size`) all sit at proseg's own defaults → none is forced onto a
  sweep axis.
- Of the six previously-hardcoded knobs, only **`cell_compactness` deviated** (hardcoded 0.03 vs
  proseg's true 0.04). Per the audit rule its config default is reset to proseg's 0.04 and the
  deviating 0.03 is covered by the sweep (which straddles 0.03 and both directions).
- `--nthreads` is a pure performance/resource knob (computed from `meta['cpus']`) → left fixed,
  not swept.

**Tiers:**

- **Tier 0 — input, not a parameter.** The prior `segmentation.zarr` (nuclear/cell label
  quality) is the biggest lever — proseg refines it, it does not create cells from nothing — and
  the transcript density/panel. Not in the sweep.
- **Tier 1 — the strongest model levers (swept):**
  - `cell_compactness`: 0.04 → **0.02, 0.03, 0.06** (smaller = more compact cells; walks through
    the old hardcoded 0.03 and a looser 0.06).
  - `nuclear_reassignment_prob`: 0.2 → **0.05, 0.5** (how freely transcripts leave the nuclear
    prior label).
  - `diffusion_probability`: 0.2 → **0.1, 0.35** (how aggressively transcripts are treated as
    diffusion noise and repositioned).
- **Tier 2 — capacity / convergence knobs:**
  - `ncomponents`: 10 → **5, 15** (mixture components for the gene-expression model; swept).
  - `samples` / `burnin_samples` (200/200): sampler convergence budget — more buys accuracy at a
    large runtime cost. At proseg default, exposed but NOT swept here (Tier-2, keeps the test
    sweep light); candidates for a dedicated quality-vs-time sweep.
- **Tier 3 — exposed but not swept (interaction with `diffusion_probability`):**
  `diffusion_sigma_far` (4.0) and `diffusion_sigma_near` (1.0) — repositioning stddevs. Left at
  default because `diffusion_probability` is the master switch already on a Tier-1 axis; useful
  ranges for a later sweep are ~2–8 (far) and ~0.5–2 (near).
- **Tier 3 — spatial-resolution knobs (exposed, not swept):** `voxel_size` (1.0) and
  `voxel_layers` (4) change the lattice resolution / 3D modelling; at proseg defaults, heavy on
  runtime/memory, left fixed for the smoke sweep.

**Total variants = 1 default + (3 + 2 + 2 + 2) = 10.**
