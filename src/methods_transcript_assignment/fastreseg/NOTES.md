# fastReseg — transcript assignment

## What this component is

fastReseg is a **transcript-assignment** method (API: `src/api/comp_method_transcript_assignment.yaml`;
inputs `raw_ist.zarr` + `segmentation.zarr` + `scrnaseq_reference.h5ad`, output
`transcript_assignments.zarr`). It corrects an initial image-based segmentation using the
spatial profile of transcripts. Wraps the R package
[`Nanostring-Biostats/FastReseg`](https://github.com/Nanostring-Biostats/FastReseg)
(DOI 10.1038/s41598-025-08733-5).

**What's unusual: it is the only multi-script, bash-orchestrated component in the repo.**
It has four resources and the *first* (`orchestrator.sh`, a `bash_script`) is the Viash
entrypoint. It chains three sub-scripts through TSV/CSV files in a temp dir:

```
orchestrator.sh
  → input.py   (SpatialData zarr  → counts.tsv, transcripts.tsv, cell_types.tsv)
  → script.R   (FastReseg::fastReseg_full_pipeline → cell_ids.csv, gene_names.csv, transcripts_out.csv)
  → output.py  (CSVs → transcript_assignments.zarr)
```

## orchestrator.sh — control flow

- `set -eo pipefail` (line 6) — **load-bearing.** Without it, a crash in any sub-step is
  swallowed (the script's last line is `echo $(date)` → exit 0) and Viash reports only the
  generic *"Required output file is missing. Expression: par.required"*, hiding the real
  Python/R traceback. This misdirection cost a whole debugging round; keep it.
- Sub-scripts are invoked as `"$meta_resources_dir/<script>"` (lines 47/60/68), **not** bare
  `python "input.py"`. The orchestrator runs from a temp copy in `$VIASH_META_TEMP_DIR`, not
  the resources dir, so bare paths work under `viash test` (cwd happens to match) but break in
  Nextflow.
- Intermediates go in `mktemp -d -p "$(pwd)"` (the writable task dir).

## input.py — step by step

1. Read `raw_ist` + `segmentation` zarrs; assign each transcript to a segmentation cell id by
   indexing the label image at transcript (y,x), handling a `global` translation if the
   segmentation was cropped.
2. Build `transcripts.tsv` (renamed cols `target`/`UMI_transID`/`UMI_cellID`) — line ~132.
   **`transcripts.compute().copy()` is load-bearing** (line 132): for a single-partition
   points object, `.compute()` returns a reference to the dask graph's backing frame, so the
   in-place `.rename()` two lines later would mutate it, and the *later* `transcripts.compute()`
   fed to `run_ssam` would yield renamed columns that no longer match the dask `_meta`
   → "Metadata mismatch". The `.copy()` isolates the rename.
3. Aggregate per-cell counts with a **local `generate_adata`** (defined at line 13), *not*
   `tx.preprocessing.generate_adata`. Copied verbatim from
   `src/methods_count_aggregation/basic_count_aggregation/script.py`: txsim's version calls
   `ad.AnnData(X, dtype=...)` (the `dtype=` kwarg was removed in anndata ≥ 0.11; image ships
   0.13) and uses unsupported per-row item assignment. Keep the two copies in sync.
4. Cell-type annotation via `tx.preprocessing.run_ssam` (line ~169) using the SC reference's
   `normalized` layer → `ct_ssam` → `cell_types.tsv`. `run_ssam` pulls in `plankton`
   (see Setup — the matplotlib pin exists for this).

## script.R — step by step

`FastReseg::fastReseg_full_pipeline(...)`. Note the branch at the end: if
`updated_transDF_list[[1]]` is NULL (**"no transcripts assigned after refinement"**), it falls
back to a `left_join` of the ssam cell types onto the input transcripts and renames
`ct_ssam → updated_celltype`, `UMI_cellID → updated_cellID`. On the small
`mouse_brain_combined` test data this fallback *is* taken and cell types come back mostly
`None_sp` — the component still produces valid output, but whether FastReseg does anything
meaningful at that data scale is unverified (see Risk points).

## output.py — CSVs → zarr

Builds a SpatialData with points `transcripts` (`x,y,z,feature_name,cell_id,updated_celltype,
transcript_id`) and a `table` whose `obs` carries `cell_id`/`cell_type` (renamed from
FastReseg's `updated_cellID`/`updated_celltype` — required by
`src/api/file_transcript_assignments.yaml`).

- **`feature_name` is cast to `category` (line 44) — load-bearing for the *downstream* stage.**
  Rebuilt from R's CSV, `feature_name` is otherwise a pyarrow-backed nullable `string` that
  survives the zarr round-trip. The next stage, `basic_count_aggregation`, builds `var_names`
  from `pd.unique(feature_name)` and `write_h5ad` then refuses the nullable-string `var`
  `_index` (*"anndata.settings.allow_write_nullable_strings is None"*). `raw_ist` transcripts
  are `category`; matching that convention fixes it. **This is not caught by fastReseg's own
  `viash test`** — that test only checks the transcript-assignment output *format*, not the
  downstream consumer.

## Setup / Docker (highest-value section)

Base `openproblems/base_python:1` (**Debian 13 trixie, amd64-only**) + merges
`setup_txsim_partial` + `setup_spatialdata_partial`. Setup step ORDER matters — Viash runs
merged partials first, local `setup:` last, so:

1. (merged) txsim: squidpy, rasterio, `theislab/txsim@dev`
2. (merged) spatialdata: `spatialdata>=0.7.3`, `anndata>=0.12.0`, `zarr>=3.0.0`
3. (local, apt) system libs + `r-base` + **FastReseg's R deps as prebuilt `r-cran-*` Debian
   binaries** (lines 55–72)
4. (local, docker run) `remotes::install_github("Nanostring-Biostats/FastReseg",
   upgrade = "never")` (line 79)
5. (local, python) **`planktonspace` + `matplotlib<3.9`** (line 86) — runs LAST so it beats
   the `matplotlib ≥ 3.9` that spatialdata/squidpy pull in.

Why each non-obvious choice:

- **`matplotlib<3.9`** — `plankton` (imported by `txsim.run_ssam`) does
  `from matplotlib.cm import get_cmap`, removed in matplotlib 3.9. Mirrors
  `methods_cell_type_annotation/ssam`, which pins the same. (Side effect: drags scanpy down
  to 1.11.)
- **`r-cran-*` binaries instead of source compilation** — the R stack (terra, igraph, sf,
  spatstat, stringi, …) compiled from source under amd64 emulation took ~1.5–2 h. Debian
  trixie ships these as prebuilt binaries; installing via apt is near-instant even emulated,
  cutting the R stage to ~15–30 min. Everything FastReseg needs is packaged **except**
  `concaveman`, `GiottoUtils`, `GiottoClass` (and FastReseg itself) — those still compile,
  but their heavy deps are already present.
- **`upgrade = "never"`** — without it `remotes` re-fetches and recompiles the apt-installed
  binary deps from CRAN, defeating the whole optimization.
- **`devtools` was dropped** — it was in the R list but is *not* a FastReseg dependency and
  isn't used at runtime; it dragged in a large subtree (roxygen2/pkgbuild/gert/…).

## Arguments

Only the standard API args (`--input_ist`, `--input_segmentation`, `--input_scrnaseq`,
`--sc_cell_type_key`, `--output`). Method-specific knobs (pixel/z-step size, distance cutoffs,
lrtest/svm cutoffs) are currently **hardcoded** in `script.R` and `um_per_pixel`/
`sc_celltype_key` in `input.py` — not yet exposed as Viash arguments (the commented
`arguments:` block in the config is a placeholder).

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: dep `methods_transcript_assignment/fastreseg`
  (line 170); appended to the `--transcript_assignment_methods` default (line 71).
- `src/workflows/run_benchmark/main.nf`: in `segm_ass_methods` (line 143), fanned out via
  `expandChannelWithParameterSets(..., "ass", ...)`.

## Risk points / gotchas

- **Two copies of `generate_adata`** (here + `basic_count_aggregation`) must stay in sync.
- **`txsim@dev` is unpinned** — its `generate_adata`/`run_ssam` contracts can drift; that's
  exactly what broke here. Consider commit-pinning.
- **`updated_celltype` in the output points** is extra (beyond the API's cell_id/feature_name);
  harmless because `basic_count_aggregation` reads only cell_id + feature_name, but non-standard.
- **Science unverified**: the "no transcripts assigned after refinement" fallback + `None_sp`
  cell types on test data mean the method's actual behaviour on real data is not yet confirmed.
- **Validated**: `viash test` passes locally (amd64/emulated); all four runtime bugs and the
  downstream `basic_count_aggregation` compatibility were reproduced and fixed in a k8s pod on
  the `build_main` image. Not yet run through a full Nextflow benchmark end-to-end.
