# bruker_cosmx loader — notes

## What this component is

Dataset loader (`datasets/loaders` namespace) that converts a **Bruker / NanoString
CosMx** export into the standardized `raw_ist` SpatialData zarr consumed by the rest of the
pipeline (API: `/src/api/file_common_ist.yaml`, merged as `--output` in `config.vsh.yaml`).

It is a thin wrapper around **`sopa.io.cosmx`** (the sopa dependency is unpinned in
`config.vsh.yaml` → currently sopa ~2.2.x). sopa does the heavy lifting: it stitches the
per-FOV morphology images and cell-label tifs into one global image/label, reads the
transcripts flat file, and builds the cells table. Our script's job is (1) massage the
on-disk directory layout into the shape sopa expects, (2) rename sopa's element keys to our
API names, (3) drop non-DAPI image channels, (4) attach dataset metadata, (5) write.

Some CosMx exports ship **no AtoMx flat files at all** (the liver "Raw Data Files" release —
see Version 4 below). For those the script first **reconstructs** the flat files sopa needs
from the raw per-FOV decoding CSVs, then feeds the same `sopa.io.cosmx` call. This is
detected automatically (no config knob): if `input_flat_files` is unset AND the raw zip has
no `*_tx_file.csv`, reconstruction mode is on.

Sibling loader `bruker_cosmx_nsclc` handles the CosMx lung-cancer layout (per-sample
subdirs, per-z-plane images) — that's a *different* export structure (see the docstring
"Version 3"), not covered here.

Datasets driven through this loader (see `scripts/create_resources/spatial/process_bruker_cosmx_nebius.sh`):
- `bruker_cosmx/bruker_mouse_brain_cosmx/rep1` — full mouse-brain hemisphere (`HalfBrain.zip`
  + separate `Half Brain simple files.zip`). **The big one** — the source of the OOM this
  NOTES documents.
- `bruker_cosmx/bruker_human_liver_cosmx` — human liver (`NormalLiverFiles.zip`). **No flat
  files exist for this dataset** (NanoString only released raw images + TileDB + Seurat), so
  it runs through the **reconstruction** path (Version 4 / see section below). `input_flat_files`
  is left unset in the run script.

## Supported input layouts

The script's module docstring (`script.py:3-80`) is the authoritative catalogue of the three
CosMx export shapes seen in the wild. The two this loader handles:

- **Flat files bundled in the raw zip** (liver): `input_flat_files` is left unset.
- **Flat files in a separate zip** (mouse brain): `input_flat_files` is set, and the CSVs are
  symlinked into `CellStatsDir` (`script.py:170-184`).

In both, sopa wants a `CellStatsDir/CellLabels/` folder, but many exports only ship the
per-FOV `CellLabels_F###.tif` inside `FOV*/` folders. `script.py:186-200` gathers those into
a `CellLabels/` folder via symlinks (see sopa issue #285).

## script.py — step by step

1. **Extract** `input_raw` (strip the single wrapper dir) → find `CellStatsDir`
   (`script.py:155-168`).
2. **Extract flat files** if provided, symlink the `*.csv` into `CellStatsDir`
   (`script.py:170-184`).
3. **Assemble `CellLabels/`** from per-FOV tifs if needed (`script.py:186-200`).
4. **Lean transcript-read shim** installed here (`script.py:202-253`) — see next section.
5. **`sopa.io.cosmx(...)`** with `cells_labels/table/polygons=True`, `read_proteins=False`
   (`script.py:261`). Stitched image/labels come back as **lazy dask** — they are *not* the
   memory peak; they're written chunked.
6. **Rename element keys** to API names (`script.py:277-288`):
   `stitched_image→morphology_mip`, `stitched_labels→cell_labels`, `points→transcripts`,
   `cells_polygons→cell_boundaries`.
7. **Add API transcript columns** (`script.py:290-292`): `cell_id = global_cell_id`,
   `feature_name = target` (dask-lazy copies; cheap now that `target` is categorical).
8. **Keep only the DNA channel** of the morphology image (`script.py:300`,
   `.sel(c=["DNA"])`) — we treat the "DNA" stain as DAPI. TODO in code: keep PolyT /
   Cellbound channels one day (currently saving/plotting fails when they're kept).
9. **Attach dataset metadata** to `table.uns` (`script.py:303-307`).
10. **`sdata.write(par["output"])`** (`script.py:315`).

## Flat-file reconstruction for raw-only exports (Version 4, liver — 2026-07-29)

`NormalLiverFiles.zip` contains only imaging + segmentation + raw decoding output; the AtoMx
"flat files" were never released for the liver dataset (verified: no `*_tx_file.csv` /
`*_fov_positions_file.csv` anywhere in the zip, and no companion "simple files" zip on S3 —
only TileDB/Seurat). `sopa.io.cosmx` therefore died at `_infer_dataset_id`
(`ValueError: Could not infer dataset_id ...`). The raw material to rebuild the flat files is
present, so the loader now reconstructs them.

**Detection** (`script.py`, `RECONSTRUCT`): on when `input_flat_files` is unset *and* the raw
zip has no `*_tx_file.csv` (peeked cheaply from the zip central directory via
`_zip_has_flat_files`). Off for mouse brain (separate flat-files zip) and any export that
ships flat files inside the raw zip.

**Extraction** switches from the denylist to an *allowlist* (`_member_wanted_reconstruct` /
`_RECON_KEEP`) that keeps only what reconstruction + sopa need: per-FOV `CellLabels_F*.tif`,
`*Cell_Stats_F*.csv`, `Morphology2D/*.TIF`, `AnalysisResults/**/*complete_code_cell_target_call_coord.csv`,
and `RunSummary/*.fovs.csv`. Footprint for the liver hemisphere ≈ **61 GB** (vs 440 GB of
Morphology3D alone). Crucially this keeps `AnalysisResults`/`RunSummary`, which the flat-file
denylist drops.

**`reconstruct_flat_files()`** writes four flat files into `CellStatsDir` (`DATA_DIR`):
- `{id}_fov_positions_file.csv` ← `RunSummary/latest.fovs.csv` (headerless; `x_mm`=col 1,
  `y_mm`=col 2, `fov`=col 6 — matches the Metrics4D `Xpos_um`/`Ypos_um`).
- `{id}_tx_file.csv` ← concat of the 304 `complete_code_cell_target_call_coord.csv`
  (`fov, cell_ID=CellId, target, x_global_px, y_global_px, z`), written FOV-by-FOV so only
  one FOV is in RAM (≈229 M transcripts total for liver).
- `{id}_metadata_file.csv` ← `Cell_Stats_F*.csv` (`fov, cell_ID, CenterX/Y_global_px, Area`).
- `{id}_exprMat_file.csv` ← per-cell gene counts crosstab'd from the assigned transcripts,
  reindexed to the Cell_Stats cell list so metadata and counts share the exact same cells/order.
- No `-polygons.csv`: the raw export has no cell boundaries → `cells_polygons=False`,
  `cell_boundaries` absent (optional in the API; downstream methods segment themselves).

**Coordinate formula** (the load-bearing part — validated end-to-end at **0.998** of assigned
transcripts landing on their own cell in sopa's stitched labels, on a 3-FOV subset run through
the real `sopa.io.cosmx`):
```
scale       = 1e3 / COSMX_PIXEL_SIZE      # 0.120280945 µm/px; == sopa's own mm→px factor
x_global_px = x_mm * scale + x_local
y_global_px = y_mm * scale + (H - 1 - y_local)   # H=4256; sopa flips each label tile along y
```
The transcript local `(x, y)` map DIRECTLY to `CellLabels[row=y, col=x]` (top-origin, no flip;
measured 0.993 raw), but sopa's stitcher does `da.flip(tile, axis=y)`, so the `(H-1-y)` term
undoes that flip. The same transform is applied to the Cell_Stats centroids so the table's
`obsm['spatial']` stays consistent with transcripts + labels. `fov_shift` is left at sopa's
inferred default (False — no polygons to infer from).

**Global-cell-id assert workaround:** sopa derives `max(cell_ID)` separately from each of tx /
metadata / exprMat and asserts they agree; a segmented cell with zero transcripts makes the
tx-file max smaller and trips the assert. In reconstruction mode we monkeypatch
`_CosMXReader._get_global_cell_id` to use one fixed max (the segmentation's true max from
Cell_Stats), keeping ids consistent across all three files.

**Not yet run end-to-end on the full 304-FOV liver dataset on Nebius** — validated on a 3-FOV
subset. Watch memory (≈229 M-row tx read, mitigated by the lean-read shim) and wall-time
(reconstruction + stitch + zarr write) against the loader's `hightime` (8 h) cap; bump to
`veryhightime` if it times out.

## The real OOM: scratch exhaustion during extraction (the actual fix, 2026-07-27)

**Symptom:** `bruker_mouse_brain_cosmx` (full hemisphere) is OOM-killed (**exit 137**) on
every retry at the *very first* step — the `log("Extract zip of raw files")` line, 58 µs in —
**before any sopa / transcript / image code runs**.

**Root cause (measured — NOT the transcript read):** the raw `HalfBrain.zip` is **208 GB
uncompressed** (3990 files), but sopa only ever reads **Morphology2D (13.6 GB) +
`FOV*/CellLabels` (0.12 GB)**. ~190 GB of the zip is `Morphology3D` (per-z-plane 3D stacks —
sopa uses the 2D composite), `AnalysisResults` (decoding CSVs), and `RunSummary` — none of
which sopa touches. The loader extracted all of it. And because `--input_raw` was a Viash
`file`, Nextflow first **staged the 187.7 GB zip** into the task's scratch, then the script
extracted ~208 GB *alongside* it. The `veryhighmem` node has only ~220 GB ephemeral disk and a
**RAM-backed scratch** (which is why a *disk* operation OOMs as *memory*): staged zip (187) +
extraction blow the 400→480 GB memory limit → 137. Retries can't help — `disk` is pinned at
200 GB in `labels_nebius.config` and memory only reaches the 480 GB cap.

**Fix (2026-07-27):**
1. `--input_raw` is now `type: string` in **both** the loader and the `process_bruker_cosmx`
   workflow config, so Nextflow no longer stages the 187 GB zip (verified: `main.nf` passes it
   through `fromState`).
2. `extract_zip` (`script.py`) streams the zip **in place** — local path via `open()`,
   `s3://` via **s3fs** (`_open_zip_source`, 32 MB blocks) — and extracts only members **not**
   under `SKIP_DIRS = {Morphology3D, AnalysisResults, RunSummary}` (`_member_wanted`).
   Footprint drops **208 GB → 14.6 GB** (916 of 3990 files); S3 transfer drops ~14×. Verified
   end-to-end against the real S3 zip (dry-run filter count + streamed a Morphology2D and a
   CellLabels tif → valid TIFFs, `strip_root` correct). New dep **s3fs** in `config.vsh.yaml`
   pypi ⇒ **the container must be rebuilt** before a run picks this up.

Denylist (not allowlist) on purpose: it drops only the three known-huge unused dirs and keeps
every small file, so it can't accidentally drop the flat CSVs the **liver** dataset ships
*inside* its raw zip.

## Secondary optimization: lean transcript read (NOT the OOM — kept anyway)

The *earlier* hypothesis (recorded here before 2026-07-27) was that the OOM was sopa's
`_CosMXReader.read_transcripts` doing one `pd.read_csv` of the whole `*_tx_file.csv`. That read
is genuinely wasteful (the mouse-brain tx file is 7 GB / ~100M rows, object-dtype strings) and
the lean shim below is worth keeping — but it was **never the OOM**: the job died at
extraction, long before transcripts are read, and the panel is only ~960 genes (so sopa's
table densification is a non-issue too). The shim is only now *reached* for the first time,
once extraction succeeds.

**Lean transcript read (`script.py`):** we replace `sopa.io.reader.cosmx.pd` with a thin proxy
(`_PandasReadCsvProxy`) that forwards every attribute to real pandas but intercepts
`read_csv`. For the transcripts file only (name contains `_tx_file.csv`), it:
- reads `target` as **categorical** (object → int codes; the single biggest saver), and
- restricts `usecols` to the columns sopa's stitched-FOV path + `_get_global_cell_id` + the
  downstream schema actually use: `fov, cell_ID, target, x_global_px, y_global_px, z`
  (dropping `cell`, `CellComp`, `x_local_px`, `y_local_px`).

It pre-reads the header (`nrows=0`) to intersect the keep-list with the real columns and only
leans the read if sopa's required columns are present — otherwise it falls through to a normal
read, so it degrades gracefully on an unexpected export. Every other sopa read (fov positions,
metadata, counts, polygons) is left untouched because their filenames don't match
`_tx_file.csv`. **All of sopa's coordinate/stitching math runs unchanged** on the leaner frame.

Why the proxy (not patching `pandas.read_csv` globally, and not a private-API reimplement):
sopa is unpinned, so we avoid touching its internals; the proxy only rebinds the module-level
`pd` name *inside sopa's cosmx module*, so no global pandas mutation and no dependence on
sopa's private `read_transcripts` signature.

Measured ~5.6× smaller (82%) transcript frame on a synthetic CosMx-shaped file; larger on the
real data. Safe to drop the vendor-specific columns because the pipeline is vendor-agnostic —
downstream methods consume the standardized API columns (`x, y, feature_name, cell_id, z`);
nothing depends on CosMx-only `cell`/`CellComp` (verified: the `df["cell"]` refs in
`methods_transcript_assignment/baysor/script_no_sopa.py` are Baysor's *own* output, not our
input).

## Gotcha: obs `cell_id` vs vendor `cell_ID` (spatialdata write, 2026-07-30)

`sdata.write()` failed for the mouse brain with `ValidationError: SpatialData contains elements
with invalid names`. spatialdata's writer forbids two keys in a table attribute that differ only
in **case**. sopa's table already carries the CosMx per-FOV-local **`cell_ID`** obs column; the
loader then adds **`cell_id`** (the global unique index — required by `file_common_ist.yaml`'s
`obs.cell_id`). `cell_ID` vs `cell_id` collide → invalid.

Fix (`script.py`, "Add info to metadata table"): rename the vendor local id to `cell_ID_local`
*before* adding `cell_id`. Note this is a **case-insensitive-key** rule, not a bad-character one
— the message's "invalid names" is a catch-all. (Aside: the mouse panel also has `/` gene names
like `Tuba1a/b/c`; those live in the var **index**, which spatialdata does **not** validate — so
they are harmless for `write` and are left untouched. Only obs/var *column* names, obsm/uns
*keys*, and element names are checked.)

## Arguments

| Argument | Required | Notes |
|----------|----------|-------|
| `--input_raw` | yes | **`type: string`** (not a staged file) — URL/path to the raw zip. Streamed in place (s3:// via s3fs) and selectively extracted; see the OOM section. |
| `--input_flat_files` | no | Second zip with the `*_<suffix>.csv` flat files, when not in the raw zip (mouse brain needs it) |
| `--segmentation_id` | default `["cell"]` | Must be exactly `["cell"]` — asserted at `script.py`; CosMx ships only the cell segmentation |
| `--dataset_*` metadata | mixed | Written into `table.uns` |

## Setup / Docker

`config.vsh.yaml`: `openproblems/base_python:1` + `__merge__` of
`/src/base/setup_spatialdata_partial.yaml` + `pypi: [sopa, s3fs]`. **sopa is unpinned** — a
load-bearing risk given the transcript shim reaches into `sopa.io.reader.cosmx` (see risk
points). **s3fs** was added 2026-07-27 for streaming the raw zip from `s3://` — adding it means
the container must be rebuilt/pushed before a run sees the extraction fix.

## Wiring

- Registered as a Nextflow module dep of the per-dataset workflow
  `src/datasets/workflows/process_bruker_cosmx/` (`main.nf` runs loader → optional
  `crop_region` → `setState`).
- Resource label: **`[veryhighmem, midcpu, hightime]`** in `config.vsh.yaml`.
- Run script: `scripts/create_resources/spatial/process_bruker_cosmx_nebius.sh` (launches via
  Seqera `tw launch` against `build/main`).
- The workflow's `crop_region` step (gated on `crop_region_min_x`) is the escape hatch to
  tile/crop the section if it still won't fit — the mouse-brain params don't set it today.

## Risk points / gotchas

- **sopa is unpinned** and the transcript fix depends on `sopa.io.reader.cosmx` existing and
  reading the tx file via the module-level `pd.read_csv`. If a future sopa refactors the
  reader (renames the module, switches to `dask`/`pyarrow`, changes the tx filename match),
  the shim silently stops leaning (falls through to a normal read → OOM returns). Consider
  pinning sopa if this bites.
- **Secondary spike, not yet fixed:** sopa's `read_tables` densifies the cell×gene matrix via
  `csr_matrix(counts.values)` (`io/reader/cosmx.py`). Negligible for a ~1k panel; large if the
  mouse-brain export is **WTx (~19k genes)**. Can't fix without reimplementing sopa's table
  read — next lever if OOM persists after the transcript fix.
- **Script-only change:** no image rebuild needed, but it must reach `origin/main` and be
  regenerated into `build/main` before a Nebius run picks it up (use the `check-component`
  skill to confirm deploy-freshness).
- **Not yet confirmed end-to-end:** the fix passes local syntax + logic tests; it has **not**
  yet been validated on a real mouse-brain run on Nebius.
