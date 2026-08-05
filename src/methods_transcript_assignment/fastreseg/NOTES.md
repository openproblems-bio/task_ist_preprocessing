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

Standard API args (`--input_ist`, `--input_segmentation`, `--input_scrnaseq`,
`--sc_cell_type_key`, `--output`) **plus four exposed FastReseg tuning knobs** (added for the
parameter sweep — see "Optimization / tuning" below):

| Arg | Type | Default | Maps to (`fastReseg_full_pipeline`) |
|-----|------|---------|-------------------------------------|
| `--molecular_distance_cutoff` | double | 2.7 | `molecular_distance_cutoff` |
| `--flagCell_lrtest_cutoff` | double | 5 | `flagCell_lrtest_cutoff` |
| `--svmClass_score_cutoff` | double | -2 | `svmClass_score_cutoff` |
| `--cutoff_spatialMerge` | double | 0.5 | `cutoff_spatialMerge` |

Each default equals FastReseg's TRUE package default (verified against the upstream function
signature), so the exposed-arg behaviour with defaults is byte-for-byte the old hardcoded
behaviour. **Still hardcoded** (not exposed): `pixel_size` (0.18) and `zstep_size` (0.8) in
`script.R` (Tier-0 dataset-physical constants); `groupTranscripts_method` ("dbscan") /
`spatialMergeCheck_method` ("leidenCut") (Tier-2/3 categoricals); the NULL auto-cutoffs
(`cellular_distance_cutoff`, `score_baseline`, `lower/higherCutoff_transNum`); and
`sc_celltype_key` ("cell_type") in `input.py`.

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

## Optimization / tuning

FastReseg corrects an image-based segmentation from the spatial profile of transcripts, so its
levers control **how aggressively transcript groups are formed, flagged as mis-segmented, split,
and re-merged**. Ranges below are grounded in the upstream `fastReseg_full_pipeline()` roxygen
docs (DOI 10.1038/s41598-025-08733-5).

**Non-default audit (important):** every value that was hardcoded in `script.R` *exactly*
matched FastReseg's package default — `pixel_size=0.18`, `zstep_size=0.8`,
`molecular_distance_cutoff=2.7`, `flagCell_lrtest_cutoff=5`, `svmClass_score_cutoff=-2`,
`groupTranscripts_method="dbscan"`, `spatialMergeCheck_method="leidenCut"`,
`cutoff_spatialMerge=0.5`, and all NULL auto-cutoffs. So **nothing was a pre-tuned deviation**;
each exposed config default is set to that same package default and the sweep walks the knob
away from it in both directions.

### Tiers

- **Tier 0 (input / not swept):** `pixel_size` (0.18) and `zstep_size` (0.8) — dataset-physical
  geometry that converts transcript pixel coords to microns; wrong values rescale *all* distance
  cutoffs at once. Left hardcoded because they are properties of the data, not quality dials.
  (Caveat worth revisiting: `input.py` feeds FastReseg coords already in the standardized global
  space, so whether 0.18 is the right multiplier here is unverified — but that is a correctness
  question, not a sweep axis.)
- **Tier 1 (swept — highest impact on output quality):** the four exposed knobs below.
- **Tier 2/3 (exposed candidates, NOT swept):** `groupTranscripts_method`
  (`dbscan`↔`delaunay`) and `spatialMergeCheck_method` (`leidenCut`↔`geometryDiff`) are
  categorical algorithm choices — worth a future 1-value-each flip, but left hardcoded to keep
  this sweep focused on the continuous error-detection/correction dials. The NULL auto-cutoffs
  (`cellular_distance_cutoff`, `score_baseline`, `lower/higherCutoff_transNum`) are
  data-adaptively computed by FastReseg; pinning them to fixed values is a deeper study.

### Swept knobs (Tier 1) — 11 variants (1 default + 10 sweep)

| Arg | Default | Sweep | Rationale |
|-----|---------|-------|-----------|
| `--molecular_distance_cutoff` | 2.7 | `[1.5, 2.0, 3.5, 4.0]` | Sharpest lever: max molecule↔molecule µm distance for grouping transcripts into candidate cells. Smaller → tighter/fragmented groups; larger → looser groups that may merge neighbours. Widest sweep. |
| `--flagCell_lrtest_cutoff` | 5 | `[3, 8]` | `-log10 p` cutoff to flag mis-segmented cells. Lower → flags more cells (aggressive re-segmentation); higher → conservative. |
| `--svmClass_score_cutoff` | -2 | `[-3, -1]` | Transcript-score boundary between high/low SVM classes; straddles the default. |
| `--cutoff_spatialMerge` | 0.5 | `[0.3, 0.7]` | Fractional (0–1) spatial constraint for accepting a group-merge; lower stricter, higher more permissive. |

Defaults omitted from each sweep list (the "default" variant already covers that point). Booleans
would get exactly one flip — none of these are boolean.

### Arg-ordering contract (CRITICAL — this component is multi-script)

The knobs are passed to `script.R` as **positional args appended AFTER the 6 file paths**, so the
order must be identical in three places or the run silently corrupts (a value lands on the wrong
parameter with no error):

```
config.vsh.yaml arguments:  --molecular_distance_cutoff, --flagCell_lrtest_cutoff, --svmClass_score_cutoff, --cutoff_spatialMerge
orchestrator.sh Rscript:    ... transcripts_out.csv  $par_molecular_distance_cutoff  $par_flagCell_lrtest_cutoff  $par_svmClass_score_cutoff  $par_cutoff_spatialMerge
script.R:                   args[7]=molecular_distance_cutoff  args[8]=flagCell_lrtest_cutoff  args[9]=svmClass_score_cutoff  args[10]=cutoff_spatialMerge
```

`script.R` reads them via a `num_arg(i, default)` helper (`as.numeric`, with a fallback to the
FastReseg default if the arg is absent/empty, so an older orchestrator degrades to stock
behaviour instead of injecting `NA`).

### NEEDS REBUILD (unlike the six sibling sweeps)

These four args did **not** exist before this change — they were hardcoded literals inside the
`fastReseg_full_pipeline(...)` call. The `build/main` container run by
`run_test_fastreseg_nebius.sh` (`--revision build/main`) will **reject** them with "unknown
option" until the component is rebuilt. Before launching the sweep: commit + push the config /
orchestrator.sh / script.R edits, run `viash ns build`, rebuild the fastreseg container on ghcr
(`build_main` tag — see the `check-component` skill), then launch. Sweep files live at
`scripts/run_benchmark/param_sweep/fastreseg_params.yaml` (+ `run_test_fastreseg_nebius.sh`).
