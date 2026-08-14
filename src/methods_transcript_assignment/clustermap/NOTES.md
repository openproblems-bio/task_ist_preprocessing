# clustermap — NOTES

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the `clustermap`
transcript-assignment method. Committed with the code.

## What this component is

- **Stage / API:** `methods_transcript_assignment` (merges
  `/src/api/comp_method_transcript_assignment.yaml`). Takes `raw_ist.zarr` (image + transcript
  points) plus a `segmentation.zarr`, and writes `transcript_assignments` — a SpatialData zarr
  with the `transcripts` points frame carrying a `cell_id` column plus a `table` whose `.obs`
  lists the cell ids.
- **Underlying tool:** **ClusterMap** (He et al. 2021, Nat Commun), installed from GitHub
  (`wanglab-broad/ClusterMap`, `pip install -e .` into `/ClusterMap`, added to `sys.path`).
  It is a *segmentation-free* method: it clusters the RNA spots themselves (density-peak
  clustering / DPC over neighbour-gene-composition + optional DAPI points) into cells. Here it
  is used as a transcript-*assignment* method — the segmentation input is required by the API
  but ClusterMap re-derives cells from the spots + the raw DAPI image.
- **GPU:** none, CPU-only. Nextflow directives `[ veryhightime, midcpu, highmem ]` (it is slow
  and memory-hungry — hence the tiled driver, see below).
- **Links:** docs/repo https://github.com/wanglab-broad/ClusterMap .

## Citation

`references.doi: 10.1038/s41467-021-26044-x` is **correct**: He, Bhatt, et al. 2021,
*"ClusterMap for multi-scale clustering analysis of spatial gene expression"*, Nat Commun 12,
5909. This is the method-specific paper (not a generic framework citation). No caveat needed.

## script.py — step by step

1. Read `input_ist` and `input_segmentation` zarrs; assert the `coordinate_system` exists on
   both the transcripts and the segmentation.
2. Transform the transcript points to the pixel `coordinate_system`, undoing any segmentation
   translation (crop handling). Multi-partition parquet indices are rebuilt to a clean single
   partition first (duplicate-index guard).
3. Extract the DAPI image: `dapi_image = np.squeeze(sdata['image']['scale0']['image'])` — the
   raw image, **squeezed to 2D (Y,X)** for these datasets, so `num_dims = 2`.
4. `format_spots_for_clustermap` — rename to ClusterMap columns
   (`spot_location_1/2/3`, `gene_name`), cast coords to int, encode gene names to 1-based ints.
5. Run ClusterMap. Because `window_size` defaults to `700` (not `None`),
   **`run_clustermap_over_chunks` is the path taken**: build a full-image model, tile the DAPI
   image into `window_size`-sized windows (`get_img`/`split`, margin = 10% of window), and for
   each tile with >=20 spots instantiate a per-tile `ClusterMap`, `preprocess()`, set
   `min_spot_per_cell`, `segmentation()`, then `stitch()` the tile back. (The un-tiled
   `run_clustermap` path only runs if `window_size is None`.)
6. Shift labels by +1 (ClusterMap uses -1 for background → 0), attach `cell_id` to the
   transcripts frame, build the `table.obs` from the unique non-zero cell ids, write the zarr.

### How each ClusterMap call receives the args (this is where behaviour lives)

Per tile the script does:
- `ClusterMap(..., xy_radius, z_radius, fast_preprocess, gauss_blur, sigma)` → constructor,
  which calls `binarize_dapi(dapi, fast_preprocess, gauss_blur, sigma)`.
- `model_tile.preprocess(dapi_grid_interval, pct_filter, LOF, contamination)`.
- `model_tile.min_spot_per_cell = min_spot_per_cell`.
- `model_tile.segmentation(cell_num_threshold, dapi_grid_interval, add_dapi, use_genedis)`.

## Arguments (exposed in config.vsh.yaml) + INERT args

Many config descriptions are `"?"`; the table below is derived by reading the script + the
upstream `ClusterMap`/`preprocessing.py` source, and records the **ClusterMap true defaults**.

| config arg | type | config default | ClusterMap true default | reaches the tool? | notes |
|-----|------|----------------|-------------------------|-------------------|-------|
| `transcripts_key` | string | `transcripts` | n/a | — | structural (points key) |
| `coordinate_system` | string | `global` | n/a | — | structural (pixel CS) |
| `window_size` | int | `700` | n/a (script-only tiling) | yes (driver) | **PERF/MEMORY** — sets tile size for `get_img`/`split`. Smaller = less RAM, too small = bad cells. Not a ClusterMap knob. **Left fixed.** |
| `use_dapi` | bool | `True` | n/a | **NO — INERT** | the string `use_dapi` **never appears** anywhere in script.py. Completely unused. |
| `xy_radius` | int | `40` | `None` (must be set) | yes (ctor + preprocess radius) | neighbourhood radius in **pixels**; sharpest lever. |
| `z_radius` | int | `0` | `None` | yes (ctor) | only relevant for 3D; data is squeezed to 2D (`num_dims=2`) so it is inert here. |
| `fast_preprocess` | bool | `False` | `False` | yes (ctor→`binarize_dapi`) | faster DAPI binarisation path (Otsu), quality/speed trade. |
| `gauss_blur` | bool | `True` | **`False`** | yes (ctor→`binarize_dapi`) | **DEVIATES** from tool default. Gaussian-blurs DAPI before binarising. |
| `sigma` | double | `1.0` | `1` | yes (ctor→`binarize_dapi`) | blur strength; only active when `gauss_blur=True`. |
| `pct_filter` | double | `0.0` | **`0.1`** | yes (preprocess) | **DEVIATES**. Fraction of lowest-density spots removed. `0.0` removes **nothing** (percentile(·,0)=min). |
| `LOF` | bool | `False` | `False` | yes (preprocess) | Local-Outlier-Factor spot removal switch. |
| `contamination` | double | `0` | **`0.1`** | yes (preprocess) | **DEVIATES**, but **only read inside `if LOF:`** → dead code at the shipped `LOF=False`. `sklearn.LocalOutlierFactor` requires `contamination ∈ (0,0.5]`, so `contamination=0` would *crash* if `LOF` were flipped True. |
| `min_spot_per_cell` | int | `5` | `5` (tutorial value; class sets no default) | yes (`erase_small_clusters`) | min transcripts to keep a cell. |
| `add_dapi` | bool | `True` | `True` | **NO — INERT (BUG)** | script passes `add_dapi=par["xy_radius"]` (a copy-paste bug), so the config value is ignored and `add_dapi` is forced to `xy_radius` (truthy → always True). |
| `use_genedis` | bool | `True` | `True` | **NO — INERT (BUG)** | same bug: script passes `use_genedis=par["xy_radius"]`, so config value ignored, forced truthy True. |
| `dapi_grid_interval` | int | `5` | `5` | yes (preprocess + segmentation) | grid spacing of DAPI points added to DPC; smaller = denser/slower. |
| `cell_num_threshold` | double | `0.1` | **`0.01`** | yes (segmentation) | **DEVIATES** (10×). "Larger → more cells." |

### Known bugs (documented, not fixed here)

- `add_dapi` / `use_genedis` are forwarded as `par["xy_radius"]` in **both** call sites
  (script.py L250-251 and L266-267). Fixing needs a config/script edit **plus a container
  rebuild**, so it is out of scope for a submittable-as-is sweep. Flagged as future work.
- `use_dapi` is a config arg that the script never reads.

## Optimization / tuning

Goal: tune ClusterMap for assignment **quality** on the standardized iST test data (images on a
~1 µm/px grid, 2D). The sweep is a **star around the shipped defaults** (one arg varied per
variant, not a grid) consumed by `run_benchmark` via `method_parameters_yaml`. Only args that
**genuinely reach ClusterMap and change the output** are swept.

**Tier 0 — the input (not a parameter).** The raw DAPI image (squeezed to 2D) and the
transcript spots ARE the input to ClusterMap; the API's `segmentation` input is ignored by the
algorithm. Biggest implicit lever, not swept.

**Tier 1 — highest impact (swept):**
- **`xy_radius`** (default 40) — spot-neighbourhood radius in pixels; the single sharpest lever
  (controls NGC neighbourhoods + density). On the ~1 µm/px grid, 40 px ≈ 40 µm — large versus
  typical cell radii. Sweep `[20, 30, 60]` to straddle tighter and looser neighbourhoods.
- **`cell_num_threshold`** (default 0.1; **tool 0.01**) — DPC cell-count threshold, "larger →
  more cells". Deviates 10× from the tool default. Sweep `[0.01, 0.05, 0.2]` — walks back to the
  tool default (fewer/larger cells) and up (more/smaller cells).
- **`dapi_grid_interval`** (default 5) — density of injected DAPI points; smaller = denser DAPI
  influence + slower. Sweep `[3, 10]`.
- **`min_spot_per_cell`** (default 5) — minimum transcripts to keep a cell; higher deletes small
  clusters. Sweep `[10, 20]`.

**Tier 2 — quality knobs / must-sweep deviations (swept):**
- **`pct_filter`** (default 0.0; **tool 0.1**) — low-density spot removal. Shipped `0.0` removes
  nothing; deviation forced onto an axis. Sweep `[0.05, 0.1]` (toward the tool default).
- **`gauss_blur`** (default True; **tool False**) — DAPI blur before binarising. Boolean
  deviation → one non-default value. Sweep `[false]` (the tool default; the True baseline is the
  default variant).

**Tier 2 — considered but NOT swept (rationale):**
- **`sigma`** (1.0) — blur strength; matches tool default and is a secondary companion to
  `gauss_blur`. Dropped to keep the variant budget in range; available for a follow-up sweep.
- **`fast_preprocess`** (False) — a speed/quality trade (approximate binarisation). Matches tool
  default; dropped for budget. Would be a `[true]` boolean flip.
- **`z_radius`** (0) — 3D-only; data is 2D (`num_dims=2`), so it is inert. Not swept.
- **`LOF` / `contamination`** — a **coupled, unsafe** pair in the star model: `contamination`
  (0 vs tool 0.1) is dead code while `LOF=False`, so sweeping it alone changes nothing (a
  duplicate of the default variant); and flipping `LOF=True` alone crashes because
  `sklearn.LocalOutlierFactor` rejects `contamination=0`. Exercising outlier removal safely
  needs a `default:`-block change (`LOF=True` + `contamination>0`) — out of scope. Not swept.
- **`window_size`** (700) — memory/perf tiling knob, not a ClusterMap parameter. Left fixed.

**Excluded because INERT (see bugs above):** `use_dapi` (never read), `add_dapi`, `use_genedis`
(both overwritten by `xy_radius` in the script → config value never reaches the tool).

**Non-default audit result:** deviations from ClusterMap true defaults were `gauss_blur`,
`pct_filter`, `cell_num_threshold` (all placed on sweep axes) and `contamination` (dead code at
`LOF=False`; excluded with rationale above). No perf-only deviation was silently swept
(`window_size` left fixed and noted).

**SUBMITTABLE AS-IS** — every swept arg (`xy_radius`, `cell_num_threshold`,
`dapi_grid_interval`, `min_spot_per_cell`, `pct_filter`, `gauss_blur`) already exists in
`config.vsh.yaml` and is correctly forwarded to ClusterMap, so the `--revision build/main`
container accepts them with **no rebuild**. Total = **14 variants** (1 default + 13 sweep).

The sweep lives in `scripts/run_benchmark/param_sweep/clustermap_params.yaml`; launch with
`scripts/run_benchmark/param_sweep/run_test_clustermap_nebius.sh` (commit + push the params
file first — it is read from GitHub at runtime).
