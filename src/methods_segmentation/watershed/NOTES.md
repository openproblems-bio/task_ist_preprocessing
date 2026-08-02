# Watershed segmentation — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how the
classic watershed pipeline works, which of its ~40 exposed knobs actually matter,
and where the sharp edges are.

Links: docs (this repo) https://github.com/openproblems-bio/task_ist_preprocessing ·
implementation https://github.com/theislab/txsim (`txsim.preprocessing.segment_watershed`).

**Citation caveat:** `references.doi` is `10.1109/34.87344` — Vincent & Soille (1991),
"Watersheds in digital spaces: an efficient algorithm based on immersion simulations",
IEEE TPAMI. That is a legitimate *foundational* watershed reference, so this is **not** a
wrong-paper mismatch like cellposev4. Minor nuance: the code uses
`skimage.segmentation.watershed`, which is the **marker-controlled / priority-flood**
variant (Meyer's flooding, plus Neubert & Protzel for the `compactness` option), not the
Vincent-Soille immersion algorithm per se. Also of historical note: the `repository` link
once (commits `d7c5ae0a0`..`5885bdac8`) pointed at `BennyStrobes/Watershed` — a completely
**unrelated genomics** method — before being corrected to `theislab/txsim` in `5885bdac8`.
The DOI is fine; leave it.

## What this component is

One of the methods at the **segmentation** stage of the iST preprocessing benchmark
(`src/methods_segmentation/`). It implements the standard API
`src/api/comp_method_segmentation.yaml`: `raw_ist.zarr` in → `segmentation.zarr` out (a 2D
cell **label map** + a stub `table`).

Unlike the deep-learning segmenters (cellpose/cellposev4/stardist), this is a **classic,
CPU-only image-processing pipeline** with **no learned model** — a chain of skimage
operations wrapped by `txsim.preprocessing.segment_watershed`. Its distinguishing feature
is an unusually **large, fully function-selectable parameter surface** (~40 args): each
pipeline stage's algorithm is chosen by a `*_func` string, and each stage's parameters are
individually exposed. Most of that surface is low-value plumbing; see "Optimization /
tuning" for the handful of knobs that move the result.

## How the watershed pipeline works (`segment_watershed`)

`img[0]` (first channel, a single 2D plane) flows through this fixed chain in txsim
(`txsim/preprocessing/_segmentation.py:382`); each stage is toggled by whether its
`*_func` string maps to a known function:

1. **Normalize** (`normalize_func` → `adjust_gamma`|`adjust_log`|`adjust_sigmoid`). Default
   `gamma` with gamma=1, gain=1 → effectively identity.
2. **Contrast** (`contrast_adjustment_func` → `equalize_adapthist`|`equalize_hist`|
   `rescale_intensity`). Default `equalize_adapthist` (CLAHE, clip_limit 0.01).
3. **Blur** (`blur_func` → `gaussian`|`median`). Default `gaussian`, `blur_sigma=1`.
4. **Threshold → binary nuclei mask** (`threshold_func` → `otsu`|`triangle`|`local_otsu`).
   Image is first cast to ubyte (`skimage.util.img_as_ubyte`; the module-level `import
   skimage` at `_segmentation.py:8` is what makes that name resolve). `local_otsu`
   (`rank.otsu`, default) is **adaptive** over a footprint (`threshold_footprint`=square of
   size `threshold_footprint_size`=50) → `nuclei = img >= local_threshold`. Global `otsu`/
   `triangle` pick one scalar threshold → `nuclei = img > threshold`.
5. **Mask post-processing** (loop over `post_processing_func_{i}`): default
   `remove_small_objects` (min_size 64), then `remove_small_holes` (area_threshold 64),
   applied to the **binary mask** before the distance transform.
6. **Distance transform** (`distance_transform_edt`) of the binary mask.
7. **Local-maxima markers** (`find_local_maxima` → `peak_local_max` with
   `local_maxima_min_distance`, then `label`). Each maximum becomes one watershed seed →
   one cell.
8. **Watershed** — `watershed(-distance, markers, mask=nuclei, watershed_line=True,
   **watershed_params)`. Floods basins from the markers within the nuclei mask.
9. **Background-intensity filter** (`filter_cells_based_on_local_background_intensity`,
   `_segmentation.py:330`) — if `bg_intensity_filter_bg_factor > 0` (default 0.3), drops
   cells whose mean intensity is below `bg_factor ×` local background and reindexes labels.

## script.py — step by step

1. **Build hyperparams** (`:34-38`) — `par.copy()`, convert the string `"None"` → Python
   `None` (viash passes unset optionals as the literal string `"None"`), drop `input`/
   `output`. Every remaining key is forwarded by name; `segment_watershed` reads what it
   recognises via `hyperparams.get(...)`.
2. **Lift compactness into `watershed_params`** (`:40-46`) — see "Optimization"; pops the
   exposed `--watershed_compactness` and nests it as `watershed_params={"compactness": …}`,
   the only dict `segment_watershed` forwards to `skimage.watershed`. Default 0.0 ==
   standard watershed (unchanged behaviour).
3. **Read input** (`:48-52`) — `sd.read_zarr`, `sdata['image']['scale0']` pulled full-res
   into memory (`.compute().to_numpy()`), and the `scale0` transform copied to reattach.
4. **Segment** (`:53`) — `tx.preprocessing.segment_watershed(image[0], hyperparameters)`.
   `image[0]` = first channel only.
5. **Post-process** (`:54-57`) — `convert_to_lower_dtype` downcasts the label array to the
   smallest uint that holds `max label`; wrap as `xarray.DataArray(dims=('y','x'))`,
   `Labels2DModel.parse` with the copied transform, store as `labels['segmentation']`.
6. **Carry the metadata table** (`:59-78`) — copies a stub `table` from
   `tables["metadata"]`; resolves `cell_id` in priority order explicit `obs["cell_id"]` →
   `instance_key` column → obs index; carries optional `region`. **This block is shared
   near-verbatim with the other segmentation components** (cellpose/cellposev4) — keep in
   sync if you touch it.
7. **Write** (`:80-83`) — `rmtree` any existing output, then `sd_output.write`.

## config.vsh.yaml — Docker setup

- Base image `openproblems/base_python:1` (plain CPU Python — no CUDA/torch; this method
  never uses a GPU).
- `setup` merges `/src/base/setup_txsim_partial.yaml` (needed — the whole method is
  `txsim.preprocessing.segment_watershed`) and `/src/base/setup_spatialdata_partial.yaml`.
  Commit `050fa15e1` reordered the installs ("Change order of installs in cellpose and
  watershed") — the txsim/spatialdata install order is load-bearing for a clean env.
- Nextflow directives: `[ hightime, midcpu, midmem ]` — 8 h / 15 CPU / 50 GB. Note the
  **hightime** budget: the classic pipeline (esp. `local_otsu` rank filter with a 50 px
  footprint and the O(image) background-intensity filter over large windows) is slow on
  big panels even though it is CPU-light per pixel.

## Arguments (the surface is large; most of it is inert)

~40 args across 8 stages. Two facts to internalise before reading the list:

- **`*_func` args are function selectors with tiny, closed mappings.** A value **not** in
  the mapping silently sets that stage's function to `None` (skips the stage) — and for
  `threshold_func`/`local_maxima_func` that means downstream names (`nuclei`,
  `distance_transform`, `local_maxima`) are never defined → the script **crashes** with a
  `NameError`. Valid values: `normalize_func` ∈ {gamma, log, sigmoid}; `contrast_
  adjustment_func` ∈ {equalize_adapthist, equalize_hist, rescale_intensity}; `blur_func` ∈
  {gaussian, median}; `threshold_func` ∈ {otsu, triangle, local_otsu}; `local_maxima_func`
  ∈ {find_local_maxima}; `post_processing_func_{i}` ∈ {remove_small_objects,
  remove_small_holes}.
- **`filter_params` drops any param a chosen function doesn't accept.** So many exposed
  args are no-ops for the default function (e.g. the `distance_transform_*` return/indices
  args, the `threshold_shift_*`/`hist`/`out`/`mask` args). Don't expect them to do
  anything unless you also switch the corresponding `*_func`.

The knobs that actually move the output (defaults in **bold**):

| Arg | Default | Effect |
|---|---|---|
| `--threshold_func` | **local_otsu** | foreground/background mask algorithm; see Tier 1 |
| `--threshold_footprint_size` | **50** | neighbourhood (px) for `local_otsu` only |
| `--blur_sigma` | **1** | gaussian pre-smoothing; higher = fewer spurious maxima |
| `--local_maxima_min_distance` | **5** | min px between watershed seeds = merge↔split dial |
| `--post_processing_min_size_1` | **64** | `remove_small_objects` min mask area (px) |
| `--post_processing_area_threshold_2` | **64** | `remove_small_holes` fill area (px) |
| `--bg_intensity_filter_bg_factor` | **0.3** | drop cells dimmer than 0.3× local bg; 0 disables |
| `--watershed_compactness` | **0.0** | NEW; compact-watershed shape regularity (Tier 3) |

## Optimization / tuning

Grounded in the txsim implementation and Xenium morphology geometry: at ~0.2125 µm/px an
~8 µm nucleus ≈ **37 px diameter, ~18 px radius, ~1000 px² area** — these anchor the
ranges below. Unlike cellposev4 the defaults here are not obviously *speed*-tuned; they
are a specific (and somewhat arbitrary) pipeline config, and the default `min_distance=5`
in particular tends to **over-split** real nuclei.

**Tier 0 — the input, not a parameter.** Only `image[0]` (a single channel) is segmented
(`script.py:53`). watershed here targets **nuclei**; whole-cell boundaries are not
recoverable from a nuclear channel alone. If a membrane/boundary stain exists, no threshold
tweak substitutes for it.

**Tier 1 — highest impact on quality (these are the swept knobs):**

- **`local_maxima_min_distance`** *(default 5)*. The dominant lever — it sets how close two
  watershed seeds may be, i.e. the merge↔split dial. 5 px is far below the nucleus radius
  (~18 px), so single nuclei get multiple seeds → over-segmentation. Sweep **[10, 15, 20]**
  (walking toward the radius) to merge fragments back into one cell.
- **`threshold_func`** *(default local_otsu)*. Decides *what is foreground* at all.
  `local_otsu` is adaptive (good under uneven illumination but slow, footprint 50 px);
  global **otsu**/**triangle** use one threshold (fast; `triangle` suits a dominant
  background peak, typical of fluorescence). Categorical — sweep the two alternatives.
- **`blur_sigma`** *(default 1)*. Heavier gaussian smoothing removes spurious distance-
  transform maxima (fewer over-split cells) but blurs small/dim nuclei. Sweep **[2, 3, 4]**
  px (stays sub-nucleus).
- **`post_processing_min_size_1`** *(default 64)*. `remove_small_objects` on the binary mask
  — the precision dial against debris/partial nuclei. Sweep **[128, 256, 512]**, all well
  below a full nucleus area (~1000 px²) so real nuclei survive.

**Tier 2 — quality knobs, exposed but not swept:**

- **`threshold_footprint_size`** *(default 50)* — only active when `threshold_func=
  local_otsu`; sets the local-Otsu neighbourhood. Couple it with any `threshold_func`
  change.
- **`bg_intensity_filter_bg_factor`** *(default 0.3)* — a real recall↔precision post-filter
  (drops dim cells); set 0 to disable, raise to prune harder. Its `window_size`/`bg_size`
  (1000/2000) are large and drive much of the runtime.
- **`post_processing_area_threshold_2`** *(default 64)* — `remove_small_holes`; minor.
- **`contrast_adjustment_clip_limit`** *(default 0.01)* — CLAHE aggressiveness; can help
  low-contrast panels but interacts unpredictably with thresholding.

**Tier 3 — high-value knob newly exposed by this component:**

- **`watershed_compactness`** *(NEW; default 0.0)*. `skimage.watershed`'s compactness was
  previously unreachable (only settable via the nested `watershed_params` dict that txsim
  forwards but the component never populated). Now exposed and lifted into that dict
  (`script.py:40-46`). 0.0 = standard watershed; small positive values (try 1e-3…1e-1)
  enable **compact watershed** (Neubert & Protzel) for rounder, more regular cell shapes —
  useful when watershed produces jagged/leaky basins. Left out of the default sweep to keep
  it at 12 variants and honour the Tier-1 focus. **Requires `viash ns build` + a container
  rebuild to take effect** (see "Risk points").

**Not worth touching for a benchmark:** the normalize/contrast function *choices*, all the
`distance_transform_*` return/indices/sampling args, the `threshold_shift_*`/`hist`/`out`/
`mask` args, `blur_mode`/`cval`/`preserve_range`/`truncate`, `*_connectivity_*`, `*_out_*`.
They are either inert under the default functions (dropped by `filter_params`) or pure
plumbing.

**Suggested quality-first order:** `local_maxima_min_distance` ↑ (fix over-splitting) →
`threshold_func` / `threshold_footprint_size` → `blur_sigma` ↑ → `post_processing_min_size_1`
↑ → `bg_intensity_filter_bg_factor` → `watershed_compactness` for shape.

## Risk points / gotchas

- **Newly exposed arg needs a rebuild.** `--watershed_compactness` was added to
  `config.vsh.yaml` and threaded in `script.py`; it has **no effect until** `viash ns build`
  regenerates the target and the Docker container is rebuilt (see the `check-component`
  skill). The rest of the sweep uses pre-existing args and works against the current image.
- **Unmapped `*_func` string = silent skip or crash.** Only sweep `threshold_func` among
  its three valid values; never point a selector at a name outside its mapping (Arguments
  §). `threshold_func`/`local_maxima_func` set to an unmapped value → `NameError` because
  `nuclei`/`local_maxima` never get defined.
- **`min_distance` too small over-splits; too large merges cells.** The default 5 errs
  toward over-segmentation — this is the first knob to raise if masks look fragmented.
- **Whole image loaded into RAM** (`script.py:51`) — full-res plane; part of why the label
  is `midmem` and the stage is `hightime`.
- **Only channel 0 is segmented** (`image[0]`) — other channels ignored.
- **Shared metadata block** (`script.py:59-78`) is duplicated across the segmentation
  components; fix together.
- **Runtime is threshold/filter-bound**, not model-bound: `local_otsu` (rank filter,
  footprint 50) and the background-intensity filter (windows 1000/2000) dominate.

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: listed as a dependency
  (`methods_segmentation/watershed`, `:162`) and included in the default segmentation
  method string `custom_segmentation:cellpose:cellposev4:binning:stardist:watershed`
  (`:65`).
- `src/workflows/run_benchmark/main.nf:101`: imported.
- Run scripts (this tuning task): `scripts/run_benchmark/run_test_watershed_local.sh` and
  `run_test_watershed_nebius.sh`, driven by the committed sweep
  `scripts/run_benchmark/watershed_params.yaml` (single source of truth; local reads it via
  `$REPO_ROOT`, Nebius via its raw GitHub URL — so it must be pushed before a Nebius
  launch). Otherwise watershed is commented out in the standard `run_test_*` scripts.

## References

- **Watershed transform (the DOI in config):** Vincent L. & Soille P. (1991), "Watersheds
  in digital spaces: an efficient algorithm based on immersion simulations", IEEE TPAMI
  13(6):583-598, DOI **10.1109/34.87344**.
- Implementation: `theislab/txsim` (`dev`), `txsim/preprocessing/_segmentation.py:382`
  `segment_watershed`.
- Underlying ops: `skimage.segmentation.watershed` (marker-controlled; `compactness` per
  Neubert & Protzel), `skimage.feature.peak_local_max`, `skimage.filters.rank.otsu`.
