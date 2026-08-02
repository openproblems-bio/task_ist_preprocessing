# Binning segmentation — implementation notes

Reference for anyone (human or Claude) working on this component. Binning is the
deliberate **"poor segmentation" baseline** of the iST preprocessing benchmark: it
does not look at the image at all — it lays a fixed square grid over the image
extent and calls each tile a "cell." The only knob is the tile (bin) edge length.

Links: docs/repo https://github.com/openproblems-bio/task_ist_preprocessing.

**Citation caveat:** the config's `references.doi` is currently
`10.1101/2023.02.13.528102` = Marco Salas et al., *"Optimizing Xenium In Situ data
utility by quality assessment and best practice analysis workflows"* (bioRxiv 2023,
Nilsson lab; theislab co-authors incl. Kuemmerle, Tiesmeyer, Luecken, Theis). That
is the **txsim-associated Xenium best-practices / benchmarking paper** that
benchmarks segmentation tools — an appropriate framework citation for a txsim
segmentation baseline (binning itself has no dedicated publication; it's a trivial
grid). This DOI is **fine as-is**. Note it was *corrected*: the introducing commit
`dd3dff390` had `10.1101/2024.11.07.622434`, which is an **unrelated** single-author
paper ("Integrative cell bin segmentation ... by Voronoi", Ming Lin, Nottingham) —
a wrong DOI, already fixed to the current one.

## What this component is

A method at the **segmentation** stage (`src/methods_segmentation/`), implementing
the standard API `src/api/comp_method_segmentation.yaml`: `raw_ist.zarr` in →
`segmentation.zarr` out (a 2D cell **label map** + a stub `table` carrying
`cell_id`/`region`).

Unlike every other segmentation method here (cellpose, cellposev4, stardist,
watershed), binning **ignores pixel content**. It only uses the image's **shape**
(and, when tuning in microns, its coordinate transform) to tile the field of view
into `bin_size × bin_size` squares. Transcripts are then assigned to whichever tile
they fall in by the downstream transcript-assignment stage, so each square becomes a
"pseudo-cell." This is intentionally a lower bound on segmentation quality — it
captures no cell morphology, splits and merges real cells arbitrarily.

The heavy lifting is a two-line txsim call, `txsim.preprocessing.segment_binning`.

## script.py — step by step

1. **Build hyperparameters** (`:61-65`) — `par.copy()`, map the literal string
   `"None"` → `None`, drop `input`/`output`. Leaves `bin_size` (and `bin_size_um`).
2. **Read input** (`:67`) — `sd.read_zarr(par["input"])`.
3. **Pull image + transform** (`:70-71`) — `sdata['image']['scale0'].image` →
   `.compute().to_numpy()` for the full-res `(c, y, x)` array, and `.transform.copy()`
   for the `scale0` coordinate transform (reattached to the output label map so it
   lands in the same coordinate space as the transcripts).
4. **Resolve bin size** (`:73-82`) — if `bin_size_um` is set, convert microns →
   pixels via `pixel_size_um(...)` (see below) and `round`; else use `bin_size`
   (raw pixels) directly.
5. **Tile** (`:84`) — `tx.preprocessing.segment_binning(image[0], bin_size)`.
   `image[0]` is the first channel (a 2D `(y, x)` plane); only its **shape** is
   used. The txsim function is literally:
   ```python
   x = np.floor(np.mgrid[0:n, 0:m][0] / bin_size)   # row // bin_size
   y = np.floor(np.mgrid[0:n, 0:m][1] / bin_size)   # col // bin_size
   bins = x*(np.ceil(m/bin_size)) + y + 1           # unique 1-based label per tile
   ```
   → a full-coverage label map: every pixel belongs to exactly one square, labels
   are `1..N`, no background (label 0). So there are **no gaps** — every transcript
   lands in some pseudo-cell.
6. **Downcast + parse** (`:85-88`) — `convert_to_lower_dtype` shrinks the label
   array to the smallest uint dtype that holds `max label`; wrap as an
   `xarray.DataArray(dims=('y','x'))`, `Labels2DModel.parse` with the copied
   transform, store as `labels['segmentation']`.
7. **Carry the metadata table** (`:89-108`) — copies a stub `table` from the input's
   `tables["metadata"]`, resolving the downstream-required `cell_id` in priority
   order: `obs["cell_id"]` → the table's `instance_key` column
   (`uns["spatialdata_attrs"]["instance_key"]`) → the obs index; optional `region`
   carried if present. **This block is duplicated in `cellpose/`,
   `cellposev4/` and other segmentation scripts** (the instance_key fallback was
   added for the Xenium WTA preview behind the Atera dataset). Keep the copies in
   sync. Note the table describes the *reference* cells from metadata, not the
   binning tiles — downstream counting/aggregation turns the label map into the
   actual cell × gene matrix.
8. **Write** (`:110-113`) — `rmtree` any existing output, then `sd_output.write`.

### `pixel_size_um` helper (`:28-49`) — added this tuning pass

Reads microns-per-pixel from the image element's coordinate transform:
`element.transform` (dict of coord_system → transformation) → pick `"global"` (fall
back to the first) → `to_affine_matrix(('y','x'),('y','x'))` → average of the
`|diagonal|` scale magnitudes. On failure it prints a warning and returns **1.0**
(the standardized grid — see below), so a micron sweep still produces correct
results on current data even if the transform can't be read.

## The bin_size unit gotcha (load-bearing)

`bin_size` is in **PIXELS of the `image` element**, not microns — the txsim
function divides pixel indices by it. The physical size of a bin therefore depends
on the image's pixel size.

The saving grace: this benchmark's `raw_ist.zarr` **images are rasterized onto a
standardized ~1 µm/pixel grid** (`target_unit_to_pixels: 1` in the loaders, e.g.
`src/datasets/loaders/vizgen_merscope/script.py:245`). Verified directly on the
mouse-brain Xenium test data: the `image` `scale0 → global` transform is
identity-scale + a translation (scale `[1,1,1]`), i.e. **1 px = 1 µm**, and s0 is
`[1, 1000, 1000]` (a 1 mm² crop). So on the standardized grid **`bin_size` in pixels
≈ bin size in microns**. (Contrast the *raw* Xenium morphology at ~0.2125 µm/px —
that resolution is gone by the time binning sees the rasterized `image`.)

`--bin_size_um` (added this pass, opt-in) makes this explicit and robust: give a
physical edge length and the script converts to pixels using the actual transform,
so bins are comparable even if a future dataset is gridded at a different resolution.

## config.vsh.yaml — setup

- Base image `openproblems/base_python:1`; merges `setup_txsim_partial.yaml`
  (`theislab/txsim@dev`, for `segment_binning`) + `setup_spatialdata_partial.yaml`
  (`spatialdata>=0.7.3`, `anndata>=0.12.0`, `zarr>=3.0.0`). No GPU, no model
  download — nothing exotic; it merges the base setups and works.
- Nextflow directives: `[ midtime, midcpu, midmem ]` (4 h / 15 CPU / 50 GB). This is
  generous — the actual work is one `np.mgrid` over the image and a write; the cost
  is dominated by reading the image and writing the label zarr. CPU-only.

## Arguments

| arg | type | default | maps to |
|---|---|---|---|
| `--bin_size` | integer | `30` | tile edge length in **pixels** → `segment_binning(img, bin_size)` |
| `--bin_size_um` | double | *(unset)* | tile edge length in **microns**; when set, overrides `--bin_size` after µm→px conversion |

`--bin_size_um` is **newly exposed this pass** — it needs `viash ns build` + a
binning **container rebuild** before a run can use it (see check-component). Default
behaviour (bin_size_um unset → use bin_size pixels) is byte-identical to before.

## Optimization / tuning

The only lever is the **bin edge length**. There is genuinely nothing else — the
underlying `segment_binning` takes exactly `(img, bin_size)`.

**Tier 0 — the input.** Binning uses only the image *extent* and pixel size, not
its content, so image channel/stain choice is irrelevant here (the opposite of the
real segmenters). What matters is the **pixel size** of the grid, which sets what a
given `bin_size` means physically — handled by `--bin_size_um` (see gotcha above).

**Tier 1 — `bin_size` / `bin_size_um` (the whole knob).** Sets the pseudo-cell
size. Grounded in the ~1 µm/px grid and brain/Xenium cell scale (nuclei ~5–8 µm,
whole cells ~10–15 µm):
- too **small** (≈ 10 µm) → many bins per cell → over-segmentation, transcripts of
  one cell split across tiles;
- ≈ **15 µm** → roughly one whole cell per bin → the best a blind grid can do;
- default **30 µm** → 2–3 cells per bin → coarse;
- too **large** (≈ 40–50 µm) → strong under-segmentation, several cells merged.
Because binning is a *baseline*, the point of the sweep is to bracket this from
sub-cellular to super-cellular, not to find a "great" value — it shows how much the
downstream benchmark score degrades as pseudo-cells drift from real cell size.

**Tier 3 — worth adding for a serious study:** none beyond the micron knob already
added. (If ever needed: non-square / hexagonal binning, or an offset/anchor for the
grid origin, would require changes to `segment_binning` upstream in txsim, not just
this component.)

**Recommended sweep** (see `scripts/run_benchmark/binning_params.yaml`): sweep
`bin_size_um` over `[10, 15, 20, 40, 50]` around the `30` default → 6 variants (a
"star" around the default, one value at a time). Light, as befits a one-knob
baseline.

## Risk points / gotchas

- **`bin_size` is pixels, not microns** — the whole "unit gotcha" section. Safe only
  because the benchmark grid is standardized to ~1 µm/px; use `--bin_size_um` to be
  explicit / dataset-robust.
- **`--bin_size_um` is not yet runtime-validated.** It's the one non-trivial new
  code path (reads the transform, converts). The µm→px math and the
  `pixel_size_um` transform read need a first real run to confirm (config parses
  clean via `viash config view`; the default path is unchanged). The 1.0-µm/px
  fallback means a bad transform read degrades to "treat µm as px," which is correct
  on current data.
- **Full coverage, no background.** Every pixel gets a positive label — there is no
  label 0 / background. Downstream aggregation therefore counts *all* transcripts
  into some pseudo-cell; there are no "unassigned" transcripts from the grid itself.
- **Whole image loaded into RAM** (`:70`, `.compute().to_numpy()`), like the other
  segmenters — the reason for `midmem` despite the trivial compute.
- **Shared metadata/cell_id block** (`:89-108`) is duplicated across segmentation
  scripts; fix together.
- **Only channel 0's shape** is used (`image[0]`); irrelevant here since content is
  ignored, but note it if the image is ever multi-channel with differing shapes.

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: dependency
  `methods_segmentation/binning` (`:160`), in the default segmentation string
  `custom_segmentation:cellpose:cellposev4:binning:stardist:watershed` (`:65`), and
  it is the worked example in the `--method_parameters_yaml` doc (`:125`).
- `src/workflows/run_benchmark/main.nf:99`: imported.
- Run scripts: enabled in `run_test_local.sh`, `run_full_local.sh`,
  `run_full_nebius.sh`, `run_mpii_nebius.sh`, `run_test_seqeracloud.sh`,
  `run_full_seqeracloud.sh`; commented in `run_test_nebius.sh` (the canonical
  params-file placement example). Dedicated sweep runners added this pass:
  `run_test_binning_local.sh` / `run_test_binning_nebius.sh` + the committed
  `binning_params.yaml`.

## References

- **Config DOI (appropriate):** Marco Salas S. et al., "Optimizing Xenium In Situ
  data utility by quality assessment and best practice analysis workflows", bioRxiv
  **10.1101/2023.02.13.528102** (2023) — the txsim-associated Xenium
  best-practices/benchmarking paper.
- txsim (`theislab/txsim@dev`) — `segment_binning` in
  `txsim/preprocessing/_segmentation.py`.
