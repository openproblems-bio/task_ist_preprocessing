# comseg — NOTES

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the `comseg`
transcript-assignment method. Committed with the code.

## What this component is

- **Stage / API:** `methods_transcript_assignment` (merges
  `/src/api/comp_method_transcript_assignment.yaml`). Takes `raw_ist.zarr` (image + transcript
  point cloud) plus a prior `segmentation.zarr`, and writes `transcript_assignments.zarr` — the
  transcripts points with a `cell_id` column plus a `table` whose `obs` lists the cell ids.
- **Underlying tool:** **ComSeg** (Defard et al.), a *point-cloud* segmentation / transcript
  clustering algorithm — it works directly on RNA positions (with the prior segmentation as a
  spatial prior) rather than on a membrane image, so it can recover arbitrary cell shapes. Run
  here **through sopa's wrapper** `sopa.segmentation.comseg`, which slices the field into
  transcript patches, runs ComSeg per patch (dask backend), and resolves the per-patch polygons
  into `comseg_boundaries`. Pure-Python, **CPU-only**.
- **GPU:** none. Nextflow directives `[ hightime, midcpu, highmem ]`; no `gpu` label.
- **Links:** ComSeg docs https://comseg.readthedocs.io/ ; repo https://github.com/fish-quant/ComSeg ;
  sopa wrapper https://github.com/gustaveroussy/sopa .

## Citation caveat (config DOI is WRONG)

`references.doi: 10.1038/s41592-020-01018-x` in `config.vsh.yaml` is **NOT ComSeg** — that DOI
is *"Cellpose: a generalist algorithm for cellular segmentation"* (Stringer et al., Nature
Methods 2021). It was almost certainly copy-pasted from another segmentation component at
integration time (PR #45). The **correct** ComSeg citation is:

> Defard, T., Laporte, H., Bertaux, M., et al. *"A point cloud segmentation framework for
> image-based spatial transcriptomics."* **Communications Biology** 7, 800 (2024).
> **DOI `10.1038/s42003-024-06480-3`** (bioRxiv preprint `10.1101/2023.12.01.569528`; PMID 38971915).

Recommended fix (not applied here to keep this a no-rebuild sweep): set
`references.doi: "10.1038/s42003-024-06480-3"` (Viash `references.doi` also accepts a list if you
want to keep a framework citation alongside). **Flagged, not edited.**

## script.py — step by step

1. Read `input_ist` and `input_segmentation` zarrs.
2. Rebuild the transcripts as a **single-partition** dask frame with a clean `RangeIndex`
   (multi-partition parquet restarts its index per partition -> duplicate labels break sopa's
   cell_id assignment and the final write).
3. `sd.to_polygons(segmentation)` -> prior `segmentation_boundaries` (drop the `label` column;
   `make_transcript_patches` recreates it).
4. Build a **fresh in-memory** `SpatialData` (image + transcripts + prior boundaries). This is
   load-bearing: running sopa on the disk-backed source zarr writes patch elements/cache back
   into the shared source dataset and corrupts it — see memory `comseg-shared-dataset-pollution`.
5. `sopa.make_image_patches` (`patch_width`, `patch_overlap`) then `sopa.make_transcript_patches`
   (`transcript_patch_width`, `prior_shapes_key="segmentation_boundaries"`,
   `write_cells_centroids=True`).
6. Assemble the ComSeg `config` dict (the six science knobs below + `dict_scale={x:1,y:1,z:1}`
   and `gene_column`) and run `sopa.segmentation.comseg(sdata, config)` under a capped dask
   `LocalCluster` (fork start method, `n_workers` clamped to `cpus-1`, `memory_limit` set to the
   full container RAM). See memory `comseg-dask-oom-fix` for why each of those is required.
7. `sopa.spatial.assign_transcript_to_cell(...)` -> `cell_id` on the transcripts; build a
   transcripts-only `SpatialData` + a `table` whose `obs` is the unique non-zero cell ids;
   `write(output)`.

## Arguments (exposed in config.vsh.yaml)

The six SCIENCE knobs go straight into the ComSeg `config` dict. Everything else is
perf/resource/structural and is left fixed.

| arg | type | config default | sopa/ComSeg default | class | notes |
|-----|------|----------------|---------------------|-------|-------|
| `--mean_cell_diameter` | double | `15.0` | **data-derived** (`2·√(area/π).mean()` of prior seg) | science | graph scale; component pins 15µm instead of deriving. |
| `--max_cell_radius` | double | `25.0` | **mean_diam·1.75** (~26µm at 15µm); docs example 50µm | science | max centroid->RNA association distance (µm). |
| `--alpha` | double | `0.5` | `0.5` | science | alphashape 0..1; 1 = convex hull, lower = tighter. |
| `--min_rna_per_cell` | integer | `5` | **20** | science | min transcripts to keep a cell. **Deviates (relaxed).** |
| `--norm_vector` | boolean | `false` | `False` | science | normalize per-cell expression vectors before clustering. |
| `--allow_disconnected_polygon` | boolean | `true` | **`False`** | science | multi-part cell polygons. **Deviates (flipped on).** |
| `--patch_width` | integer | `2000` | n/a | perf | image-patch width (make_image_patches). Memory tile. |
| `--patch_overlap` | integer | `50` | n/a | perf | image-patch overlap. |
| `--transcript_patch_width` | integer | `200` | n/a | perf | transcript-patch width; shrinking (200->100) ~quarters per-patch graph RAM. |
| `--n_workers` | integer | `4` | n/a | resource | dask worker cap (clamped to cpus-1). OOM safeguard. |
| `--worker_memory_limit` | string | `"0"` | n/a | resource | "0" -> full container RAM (neutralises the nanny). |
| `--gene_column` | string | `feature_name` | n/a | structural | transcript gene-name column. |
| `--transcripts_key` / `--coordinate_system` | string | `transcripts` / `global` | n/a | structural | element / coord-system selectors. |

## Optimization / tuning

**Verified upstream defaults** — from sopa's `_get_default_config` (the wrapper this component
runs through; `sopa/segmentation/methods/_comseg.py`): `mean_cell_diameter` = derived from the
prior segmentation areas, `max_cell_radius` = `mean_cell_diameter · 1.75`, `norm_vector = False`,
`alpha = 0.5`, `allow_disconnected_polygon = False`, `min_rna_per_cell = 20`. The ComSeg docs'
worked example uses `max_cell_radius = 50`, `min_rna_per_cell = 5`, `alpha = 0.5`. ComSeg's
deeper clustering knobs (co-expression graph `n_neighbors` ~40, `k_nearest_neighbors` ~10,
leiden `resolution` = 1) are **not** set by this component -> ComSeg's own defaults apply.

**Non-default audit result:** two exposed defaults deviate from sopa's tool default and neither
is a perf/resource knob, so per the phase-3 rule each is FORCED onto a sweep axis:
- `allow_disconnected_polygon`: component **true** vs sopa **False** -> sweep the flip.
- `min_rna_per_cell`: component **5** vs sopa **20** -> sweep walks 5 back toward 20.
`mean_cell_diameter`, `max_cell_radius`, `alpha`, `norm_vector` sit at/near sopa's defaults but
are the algorithm's core science levers, so they are swept too. The perf/resource/structural
knobs (`patch_width`, `patch_overlap`, `transcript_patch_width`, `n_workers`,
`worker_memory_limit`, `gene_column`, `transcripts_key`, `coordinate_system`) are left FIXED and
never enter the sweep — `transcript_patch_width` is *also* a mild quality knob (patch
boundaries), but it is primarily the OOM memory dial (memory `comseg-dask-oom-fix`), so it stays
in `default:`.

Tiers:

- **Tier 0 — input, not a parameter.** The **prior segmentation** fed in as `prior_shapes`
  (custom_segmentation upstream) is the single biggest lever — ComSeg uses it as its spatial
  prior and to derive centroids; its quality dominates the result. Also the transcript density /
  panel size and `dict_scale` (coords already in µm). Called out; not swept.
- **Tier 1 — highest impact on what gets detected (swept):**
  - `max_cell_radius`: 25 -> **20, 37.5, 50** (µm) — how far a cell reaches to claim RNA.
  - `mean_cell_diameter`: 15 -> **10, 20** (µm) — the neighborhood graph scale.
  - `alpha`: 0.5 -> **0.25, 0.75, 1.0** — polygon tightness (1 = convex hull).
  - `min_rna_per_cell`: 5 -> **10, 20** — which cells survive (walk toward sopa default).
- **Tier 2 — refinement / clustering behaviour (swept, one flip each):**
  - `norm_vector`: false -> **true** — normalize expression vectors before co-expression graph.
  - `allow_disconnected_polygon`: true -> **false** — forbid multi-part cell polygons.
- **Tier 3 — high-value ComSeg knobs NOT surfaced by the component (future work; would need a
  config arg + `viash ns build` + container rebuild — NOT in the submittable sweep):** the
  co-expression / clustering parameters ComSeg exposes but this component never plumbs into the
  `config` dict: `n_neighbors` (co-expression graph, ~40), `k_nearest_neighbors` (~10),
  `clustering_method` (leiden/louvain/kmeans), `resolution` (leiden, ~1), and
  `min_nb_rna_patch`. These change how transcripts are grouped into cells and are arguably the
  next tier worth exposing after the geometry knobs. Exposing any needs: add to
  `config.vsh.yaml arguments:`, thread `par[...]` into the ComSeg `config` dict in `script.py`,
  then `viash ns build` + container rebuild (see `check-component`).

Total variants = 1 default + (2 + 3 + 3 + 2 + 1 + 1) = **13**.

**SUBMITTABLE AS-IS:** every swept arg is already exposed in `config.vsh.yaml`, so the
`--revision build/main` Nebius run accepts them with NO rebuild.

## Wiring

- Registered in the run_benchmark workflow: dep `methods_transcript_assignment/comseg`
  (`src/workflows/run_benchmark/config.vsh.yaml` L167); present in the
  `--transcript_assignment_methods` default string L71
  (`basic_transcript_assignment:baysor:clustermap:pciseq:comseg:proseg:segger:fastreseg`).
- Sweep: `scripts/run_benchmark/param_sweep/comseg_params.yaml` (committed params) +
  `scripts/run_benchmark/param_sweep/run_test_comseg_nebius.sh` (Nebius launcher, standard CPU
  compute env `5hfmdCBxMRd4nHZaJKYEQZ`, labels `task_ist_preprocessing,test,comseg` — no `gpu`).
  The launch keeps `basic_transcript_assignment` enabled alongside `comseg` so the sweep is
  scored against the baseline assignment; every other stage stays on its single default.

## Risk points / gotchas

- **Config DOI is the Cellpose paper, not ComSeg** — see the Citation caveat above.
- **Fresh in-memory SpatialData is load-bearing** — running sopa on the disk-backed source zarr
  corrupts the shared dataset (memory `comseg-shared-dataset-pollution`).
- **dask worker cap + fork + full-RAM memory_limit are load-bearing** — without them comseg is
  OOM-killed (exit 137) or crashes on spawn re-import (memory `comseg-dask-oom-fix`). If a sweep
  variant still OOMs, lower `--n_workers` or shrink `--transcript_patch_width` (both left in
  `default:`, not swept).
- **Submittable-now constraint:** the Nebius run uses `--revision build/main`, so the sweep can
  only vary args already in the build/main container (the six science knobs). Any newly-exposed
  Tier-3 arg is silently ignored until rebuilt.
