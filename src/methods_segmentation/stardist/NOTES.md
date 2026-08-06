# StarDist2D segmentation — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how it
works, why the Docker/version pins are the way they are, what the tunable knobs are,
and where it can break.

Links: docs/repo https://github.com/stardist/stardist ·
DOI 10.48550/arXiv.1806.03535.

**Citation caveat (minor — the DOI is correct):** `references.doi` is
`10.48550/arXiv.1806.03535` = Schmidt, Weigert, Broaddus & Myers, *"Cell Detection
with Star-convex Polygons"*, MICCAI 2018. That **is** the StarDist**2D** method paper,
and this component runs `StarDist2D`, so — unlike `cellposev4` — the DOI matches the
model actually used. Only nuance: the pretrained weights (`2D_versatile_fluo`) were
released later with the library, and there is a separate 3D follow-up (Weigert et al.,
*"Star-convex Polyhedra…"*, WACV 2020) that does **not** apply here. No fix needed.

## What this component is

One of the methods at the **segmentation** stage of the iST preprocessing benchmark
(`src/methods_segmentation/`). It implements the standard API
`src/api/comp_method_segmentation.yaml`: `raw_ist.zarr` in → `segmentation.zarr` out
(a 2D cell **label map** + a stub `table`).

It runs **StarDist2D** — a CNN that, for every pixel, predicts a **star-convex
polygon** (a set of radial distances to the object boundary) plus an object
probability, then reconstructs instances by non-maximum suppression over those
polygons. Because objects are represented as star-convex polygons centered on each
nucleus, it excels at **round/blob-like nuclei in crowded fluorescence images** —
the DAPI-like nuclear morphology channel of iST data. It loads a **pretrained** model
(`StarDist2D.from_pretrained`), so there is no training step here.

Unlike Cellpose-SAM (up to 3 channels), StarDist2D is a **single-channel nuclear
detector**: the script feeds it `image[0]` only (see Tier 0 below).

## script.py — step by step

1. **numpy-alias shim** (`:5-15`) — restores the deprecated `np.bool`/`np.int`/…
   aliases that numpy≥1.24 removed but stardist/csbdeep still reference. Load-bearing
   (see Setup); must run **before** `import stardist`.
2. **Read input** (`:52-54`) — `sd.read_zarr(par["input"])`, then
   `sdata['image']['scale0'].image.compute().to_numpy()` pulls the full-res image
   into RAM as a `(c, y, x)` array; the `scale0` transform is copied to reattach to
   the output. (NB: the *original* 2025-07 version read `morphology_mip`; the current
   raw_ist contract exposes the morphology plane as `image` — do not revert.)
3. **Load model** (`:59`) — `StarDist2D.from_pretrained(par['model'])`. Default
   `2D_versatile_fluo`. This also loads the model's optimized `thresholds.json`
   (prob=0.479071, nms=0.3 for that model) as the fallback thresholds.
4. **Percentile normalizer** (`:64-77`) — a csbdeep `Normalizer` subclass that min-max
   scales the image to its **1st / 99.8th percentiles** (`normalize_mi_ma`), the
   recommended StarDist preprocessing but with fixed percentile bounds. `block_size`
   and `context` are derived from the image width so a large panel is processed in
   tiles; **`min_overlap` is derived from object size, not block size** (see step 6 and
   the min_overlap gotcha below).
5. **Build eval-params** (`:85-89`) — collects the newly exposed tunables
   (`prob_thresh, nms_thresh, scale`) from `par`, **dropping any that are `None`**.
   A dropped key ⇒ `predict_instances` uses the model's own optimized value (so the
   no-args call is byte-for-byte the pre-tuning behaviour). Mirrors the cellposev4
   eval-params pattern.
6. **Segment** (`:92-131`) — `model.predict_instances_big(image[0], axes='YX',
   block_size=…, min_overlap=…, context=…, normalizer=…, **eval_params)`.
   `predict_instances_big` splits the image into `block_size` blocks, calls
   `predict_instances` on each (forwarding `**eval_params` unchanged — it only
   overrides `axes/overlap_label/return_labels/return_predict`), and reassembles the
   labels into global coordinates. `image[0]` = first channel → a single 2D plane.
   **The stitching invariant is that every predicted object is smaller than
   `min_overlap`** (an object bigger than the overlap can span a block seam and can't be
   uniquely assigned → `RuntimeError: ...violates the assumption of being smaller than
   'min_overlap'`). So `min_overlap` is set from an **object-size** bound
   (`max_object_diameter`, default **192 px**), *not* from `block_size` (the old
   `block_size // 5.5` shrank it to 64 px on small panels while real blobs reached
   ~110 px → crash). `block_size` is then grown if needed to satisfy
   `min_overlap + 2*context < block_size`, and the call is wrapped in a **retry that
   doubles `min_overlap` on that specific error** so a rare oversized blob self-heals
   instead of failing the run.
7. **Post-process** (`:100-104`) — `convert_to_lower_dtype` downcasts the label array
   to the smallest uint that holds `max label`; wrap as an `xarray.DataArray`,
   `Labels2DModel.parse` with the copied transform, store as
   `sd_output.labels['segmentation']`.
8. **Carry the metadata table** (`:106-125`) — copies a stub `table` from the input's
   `tables["metadata"]`; resolves a required `cell_id` in priority order (explicit
   `obs["cell_id"]` → the table's `instance_key` column → the obs index) and carries
   an optional `region`. **This block is duplicated verbatim in `cellpose/`,
   `cellposev4/` — keep them in sync.**
9. **Write** (`:127-131`) — `rmtree` any existing output, then `sd_output.write`.

Note: the output `table` describes the *reference* cells from metadata, not the newly
predicted masks — segmentation output here is the **label image**; downstream
aggregation turns masks into a cell × gene table.

## config.vsh.yaml — the Docker setup (the highest-value section)

Base image `openproblems/base_tensorflow_nvidia:1` (GPU-capable TF). The pip pins are
**fought-over and load-bearing** — see the history before changing them:

- **`tensorflow==2.18.0`** + **`tf_keras==2.18.0`** + **`numpy>=2.0.0,<2.2.0`** +
  **`scipy<1.15.0`**. The knot: `anndata`/`spatialdata` reference `np.bool` (re-added
  in numpy **2.0**), so `numpy<2.0` makes `import anndata` crash; but TF **2.17** caps
  `numpy<2.0`. Resolution = bump TF to **2.18** (allows numpy≥2.0), and because
  csbdeep/stardist use the **Keras-2 API via `tf_keras`**, `tf_keras` must match the TF
  version. `scipy<1.15.0` is pinned for a stardist/csbdeep compatibility break.
- Even with numpy≥2.0, stardist/csbdeep still reference the **removed** scalar aliases
  (`np.bool`, `np.long`, …) at import → hence the runtime shim at the top of the script
  (`:5-15`). Both halves (the TF/numpy pin **and** the shim) are required together.
- History of this pin (git): `06d998545` add (numpy<2.0, tf loose) → `5058466a7`
  add `scipy<1.15.0` → `aa12d73fc` drop the numpy/scipy pins → `15a214eb3` (#167) the
  current tf 2.18 / tf_keras 2.18 / numpy 2.0-2.2 / scipy pin + the alias shim.
  Don't "simplify" these back.
- Dev note (in-config): on **macOS** the TF install fails; develop in a conda env
  (`conda install -c conda-forge tensorflow`), test the Docker via gh-actions.
- Nextflow directives: `[ hightime, midcpu, highmem, gpu ]` — 8 h / 15 CPU / 100 GB /
  GPU. TF uses the GPU when present; on CPU it still runs (much lighter than
  cellposev4's SAM ViT) but slower.

## Arguments

| arg | type | default | maps to |
|-----|------|---------|---------|
| `--model` | string | `2D_versatile_fluo` | `StarDist2D.from_pretrained(model)` |
| `--prob_thresh` | double | *(unset → model 0.479071)* | `predict_instances(prob_thresh=)` |
| `--nms_thresh` | double | *(unset → model 0.3)* | `predict_instances(nms_thresh=)` |
| `--scale` | double | *(unset → no rescale)* | `predict_instances(scale=)` |
| `--max_object_diameter` | integer | *(unset → 192 px)* | `predict_instances_big(min_overlap=)` |

`prob_thresh`, `nms_thresh`, `scale` were **added during tuning** (see below). They are
**optional with no static default on purpose**: leaving them unset makes StarDist use
each model's own optimized `thresholds.json` — hard-coding e.g. `0.479071` would be
wrong for a different `--model` whose optimized thresholds differ. **They need
`viash ns build` + a container rebuild to take effect** (see `check-component`).

Not exposed:
- `--n_tiles` — a pure GPU-memory tiling knob. `predict_instances_big` **already** tiles
  the image via `block_size` (≤4096 px), so `n_tiles` would only sub-tile each block for
  the forward pass; it has **no effect on segmentation quality**, only on peak memory.
  Deliberately left out of the sweep (would just pad it).
- The normalization percentiles (`1, 99.8`) are hardcoded in the script (`:76`); could
  be exposed if results look intensity-sensitive, but not a priority (see Tier 2).

## Optimization / tuning

Two things to understand before tuning:
1. Unlike `cellposev4`, the shipped defaults here are **not** aggressively speed-tuned —
   `predict_instances` runs with each model's paper-optimized thresholds. So "optimize
   for quality" here is mostly **fitting StarDist to iST nuclei** (object size + the
   recall/precision thresholds), not walking back throttled knobs.
2. Ranges are grounded in the StarDist paper/docs and the pixel geometry of the data.
   The pretrained thresholds quoted (prob 0.479071 / nms 0.3) are those of
   `2D_versatile_fluo`.

**Tier 0 — the input, not a parameter.** The script segments **`image[0]` — a single
channel only** (`:93`). StarDist2D is a single-channel **nuclear** detector (no
membrane/multi-channel mode like Cellpose), so the lever here is **which stain is fed**
and its **normalization** — feeding the crisp nuclear (DAPI-like) plane and getting the
1/99.8-percentile normalization right matters more than any threshold tweak. Its output
is **nuclei**, not whole cells; whole-cell boundaries come from downstream expansion.

**Tier 1 — highest impact on quality (all newly exposed)**

- **`prob_thresh`** *(now exposed; default = model 0.479071, range ~0..1)*. Pure
  recall↔precision dial and the most direct knob. **Lower** (e.g. 0.3–0.4) → recover
  more / dimmer nuclei (higher recall); **raise** (e.g. 0.6–0.7) → keep only confident
  detections (higher precision). In iST you usually want to capture all nuclei, so this
  is the first axis to sweep. No speed cost.
- **`scale`** *(now exposed; default = none/1.0)*. StarDist's analogue of Cellpose's
  `diameter`: rescales the image so nuclei land near the model's trained size.
  `2D_versatile_fluo` was trained on a DSB2018 subset (varied, roughly ~20–40 px nuclei);
  Xenium morphology (~0.2125 µm/px) makes an ~8 µm nucleus ≈ **38 px**, so a mild
  down/upscale can matter. `>1` upscales (small nuclei appear larger), `<1` downscales.
  Sweep both sides (≈ 0.5 … 2.0).
- **`model`** *(already exposed; default `2D_versatile_fluo`)*. The other **fluorescence**
  pretrained model is `2D_paper_dsb2018`. `2D_versatile_he` is a **brightfield H&E RGB**
  model — wrong for single-channel fluorescence — so it is NOT swept.

**Tier 1/2 — merge control**

- **`nms_thresh`** *(now exposed; default = model 0.3, range ~0..1)*. Non-max-suppression
  IoU for overlapping polygon candidates. **Lower** → more aggressive suppression
  (touching nuclei more likely merged/dropped); **raise** → keep more overlapping
  detections in crowded tissue. Secondary to `prob_thresh`/`scale`; sweep a few values
  (≈ 0.2 … 0.5).

**Tier 3 — not exposed, worth adding only for a deeper sweep**

- **Normalization percentiles** (`pmin, pmax`, currently `1, 99.8`). For panels with
  uneven illumination / low contrast the percentile bounds change what counts as
  background; expose `pmin`/`pmax` if results look intensity-sensitive.
- **`n_tiles`** — memory only, not quality (see Arguments). Not worth a sweep axis.

**Suggested quality-first sweep order:** `prob_thresh` (0.3…0.7) → `scale`
(0.5…2.0, from pixel size) → `nms_thresh` (0.2…0.5) → `model` → `2D_paper_dsb2018`.
That is exactly the sweep encoded in `scripts/run_benchmark/stardist_params.yaml`
(**13 variants** = 1 default + 4 + 4 + 3 + 1).

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: listed as a dependency
  (`methods_segmentation/stardist`, `:161`) and included in the default method string
  `custom_segmentation:cellpose:cellposev4:binning:stardist:watershed` (`:65`).
- `src/workflows/run_benchmark/main.nf:100`: imported alongside the other seg methods.
- Sweep run scripts (this tuning): `scripts/run_benchmark/run_test_stardist_local.sh`
  and `run_test_stardist_nebius.sh`, sharing the committed
  `scripts/run_benchmark/stardist_params.yaml` (local reads it directly; Nebius reads
  it from the GitHub raw URL — so it must be **pushed** before launching).
- Already enabled in the standing run configs: `run_gpu_nebius.sh`, `run_full_nebius.sh`,
  `run_mpii_nebius.sh`, `run_test_seqeracloud.sh`, `run_full_seqeracloud.sh`; commented
  out in `run_test_local.sh` / `run_test_nebius.sh` (GPU/time cost for the quick path).

## Risk points / gotchas

- **The TF/tf_keras/numpy/scipy pin quartet + the alias shim are a matched set.** Change
  one and imports break (see Setup). This is the single most fragile thing here.
- **New args need a rebuild.** `prob_thresh`/`nms_thresh`/`scale` only exist after
  `viash ns build` + a container rebuild; a stale image silently ignores them. The sweep
  has **not yet been run end-to-end** with the new args — validated only by
  `viash config view` + a script `ast.parse`.
- **`scale` through `predict_instances_big`.** It is forwarded per-block via `**kwargs`;
  block geometry (`block_size`/`context`) stays in original-pixel units, so extreme
  `scale` values interact with the block overlap. Sanity-check masks at the block seams
  if using large `scale`.
- **Only channel 0 is segmented** (`image[0]`). Fine for single-channel iST morphology;
  StarDist2D has no multi-channel mode anyway.
- **Whole image loaded into RAM** (`:53`) — full-res plane; big panels are why the label
  is `highmem`. Tiling beyond that is `predict_instances_big`'s `block_size`.
- **Shared metadata block** (`:106-125`) is duplicated across the cellpose components;
  fix all together.
- **StarDist detects nuclei, not whole cells.** Expect nuclear-sized masks; downstream
  volume/expansion stages handle cell extent.

## References

- **StarDist2D (the model this component runs, = the config DOI):** Schmidt U., Weigert
  M., Broaddus C., Myers G., *"Cell Detection with Star-convex Polygons"*, MICCAI 2018,
  arXiv **1806.03535** (DOI 10.48550/arXiv.1806.03535).
- 3D follow-up (not used here): Weigert M., Schmidt U., et al., *"Star-convex Polyhedra
  for 3D Object Detection and Segmentation in Microscopy"*, WACV 2020.
- Docs/repo: https://github.com/stardist/stardist (see `examples/other2D/
  predict_big_data.ipynb`, the source of the `predict_instances_big` + custom-normalizer
  pattern used here).
- `2D_versatile_fluo` optimized thresholds (`thresholds.json`): prob 0.479071, nms 0.3;
  trained on a subset of the DSB 2018 nuclei-segmentation dataset.
