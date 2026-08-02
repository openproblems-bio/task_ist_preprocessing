# Cellpose 4 (Cellpose-SAM) segmentation — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how it
works, why it's a *separate* component from `cellpose`, what Cellpose 4 actually
is, and where the settings come from.

Links: docs https://cellpose.readthedocs.io/en/latest/ · repo
https://github.com/MouseLand/cellpose.

**Citation caveat:** the config's `references.doi` is `10.1038/s41592-020-01018-x`,
which is the **original Cellpose paper** (Stringer et al. 2021, Nat. Methods) — the
flow-dynamics CNN, *not* the model this component runs. The model actually used
(Cellpose-SAM / `cpsam`) is described in a **different** paper (Pachitariu,
Rariden & Stringer, bioRxiv `10.1101/2025.04.28.651001`, 2025). The config DOI is
the general Cellpose reference, not the v4-specific one; treat it accordingly.

## What this component is

One of the methods at the **segmentation** stage of the iST preprocessing
benchmark (`src/methods_segmentation/`). It implements the standard API
`src/api/comp_method_segmentation.yaml`: `raw_ist.zarr` in → `segmentation.zarr`
out (a 2D cell **label map** + a stub `table`).

It runs **Cellpose 4** (a.k.a. **Cellpose-SAM**, the `cpsam` model) directly via
`cellpose.models.CellposeModel` + `model.eval`. Adapted from
`openproblems-bio/task_spatial_segmentation` (the sibling OpenProblems task) —
the two scripts are nearly identical.

### Why it exists separately from `cellpose`

There are **two** Cellpose components in this repo, and this split is deliberate:

| | `cellpose/` | `cellposev4/` (this one) |
|---|---|---|
| cellpose version | `cellpose<4.0.0` (pinned) | `cellpose>=4.0.0` |
| model | `cyto` (CNN) via `Cellpose` | `cpsam` (ViT/SAM) via `CellposeModel` |
| call path | `txsim.preprocessing.segment_cellpose(...)` wrapper | direct `model.eval(...)` |
| extra deps | needs `txsim` | no `txsim` (leaner image) |
| arg surface | ~17 args (channels, model_type, invert, do_3D, …) | 5 args (diameter, flow_threshold, niter, min_size, resample) |

History (git):
- `35061648c` "Support cellpose v4" — first attempt to make the *original*
  `cellpose` component v4-compatible (dropped `interp`, moved `resample`, noted
  `diameter` "should be None with v4").
- `8015045ab` "Use cellpose v3 instead of v4 to stay in test time limits" —
  reverted: pinned `cellpose<4.0.0`, because the v4 SAM/ViT model was **too slow**
  to finish inside the CI/`viash test` time budget on a lean/CPU runner.
- `041ab3e23` "Cellposev4 (#166)" — instead of forcing one component to straddle
  both, added **this** dedicated `cellposev4` component with a speed-tuned default
  param set (see Arguments) and registered it in the benchmark.

So `cellpose` stays as the fast-enough v3 CNN baseline; `cellposev4` is the
heavier, more-generalizing transformer model run on GPU.

## What Cellpose 4 / Cellpose-SAM is (background)

Cellpose predicts, per pixel, a **flow vector** pointing toward the center of the
cell that pixel belongs to, plus a **cell probability**. At inference each pixel
"follows the flow" for a number of iterations to a fixed point; pixels converging
to the same point become one **mask**. (This flow-dynamics core is unchanged from
Cellpose 1–3.)

**What's new in v4 (Cellpose-SAM, `cpsam`)** — from the paper itself (Pachitariu,
Rariden & Stringer, "Cellpose-SAM: superhuman generalization for cellular
segmentation", bioRxiv `10.1101/2025.04.28.651001`, v1 2025-05-01, corresponding
author C. Stringer, HHMI Janelia; CC-BY-NC). Verified against the bioRxiv abstract:
- They **adapt the pretrained transformer backbone of a foundation model (SAM)
  into the Cellpose framework** — i.e. keep Cellpose's flow-dynamics
  mask-reconstruction, swap the U-Net CNN backbone for SAM's ViT. The flow →
  fixed-point → mask machinery is unchanged from Cellpose 1–3.
- **Headline claim:** the model "substantially outperforms inter-human agreement
  and approaches the human-consensus bound" — the paper frames prior methods as
  matching inter-human agreement, and argues a human *consensus* could roughly
  halve error rates; Cellpose-SAM pushes toward that. Hence "superhuman."
- Generalization robustness was **explicitly trained in**, to: **channel
  shuffling** (order-invariance), **cell size**, **shot noise**, **downsampling**,
  and **isotropic + anisotropic blur**. This is why in practice you don't need to
  set `diameter` or pick channels.
- Positioned as a **foundation model** that drops into the existing Cellpose
  ecosystem: finetuning, human-in-the-loop training, image restoration, and 3D
  segmentation.

Consequences in the library API (from the docs / release notes, consistent with
the paper's channel- and size-invariance):
- **Trained around a mean object diameter of 30 px** (range ~7.5–120 px) and
  largely **size-invariant** → `diameter` is optional.
- **`channels` is gone** — uses the first ≤3 channels; no channel/model-type
  selection.
- **API consolidated:** the old `models.Cellpose` class (which bundled a
  `SizeModel` for auto-diameter) and `models.SizeModel` were **removed**; only
  `models.CellposeModel` remains, and `CellposeModel()` loads `cpsam` by default.
  This is why the script uses `CellposeModel`, not `Cellpose`.
- Newer 4.2+ point releases add `cpsam_v2` / DINOv3 (`cpdino`) variants, but this
  component just pins `cellpose>=4.0.0` and takes whatever default ships.

## script.py — step by step

1. **Device pick** (`:11-13`) — `torch.device('cuda' if available else 'cpu')`,
   printed. Runs on CPU if no GPU, but see the speed caveat below.
2. **Read input** (`:46-49`) — `sd.read_zarr(par["input"])`, then
   `sdata['image']['scale0'].image.compute().to_numpy()` pulls the full-res image
   into memory, and the `scale0` transform is copied to reattach to the output.
3. **Init model** (`:51-52`) — `CellposeModel(gpu=torch.cuda.is_available())`.
   No `model_type`/`pretrained_model` given → loads the default `cpsam` weights.
4. **Build eval params** (`:54`) — collects the 6 tunables
   (`diameter, flow_threshold, cellprob_threshold, niter, min_size, resample`) from
   `par`, dropping any that are `None`. Note `min_size: -1`, `flow_threshold: 0.0`
   and `cellprob_threshold: 0.0` are **passed through** (they're not `None`) — see
   Arguments for what those values mean.
5. **Segment** (`:56`) — `model.eval(image[0], progress=True, **eval_params)`.
   `image[0]` = first channel of a `(c, y, x)` array → a single 2D plane. Returns
   `(masks, flows, styles)`; only `masks` is kept.
6. **Post-process** (`:58-65`) — `convert_to_lower_dtype` downcasts the label
   array to the smallest uint that holds `max label` (uint8/16/32/64) to shrink the
   zarr; wrap as an `xarray.DataArray(dims=('y','x'))`, `Labels2DModel.parse` with
   the copied transform, store as `sd_output.labels['segmentation']`.
7. **Carry the metadata table** (`:67-86`) — copies a stub `table` from the input's
   `tables["metadata"]`. Downstream needs a `cell_id` column; it's resolved in
   priority order: explicit `obs["cell_id"]` → the table's `instance_key` column
   (`uns["spatialdata_attrs"]["instance_key"]`) → the obs index. Optional `region`
   column is carried if present. **This exact block is shared verbatim with
   `cellpose/script.py`** — the instance_key fallback was added for exports (e.g.
   the Xenium WTA preview behind the Atera dataset) that lack an explicit
   `cell_id`. Keep the two copies in sync if you touch it.
8. **Write** (`:88-91`) — `rmtree` any existing output, then `sd_output.write`.

Note: the output `table` describes the *reference* cells from metadata, not the
newly predicted masks — segmentation output here is the **label image**; the
counting/aggregation stage downstream is what turns masks into a cell × gene table.

## config.vsh.yaml — the Docker setup

- Base image `openproblems/base_pytorch_nvidia:1` (CUDA/torch) — GPU-capable.
- `setup`:
  - `pypi: cellpose>=4.0.0` — the version pin that makes this "v4".
  - `script: from cellpose.models import CellposeModel; model = CellposeModel()` —
    a **build-time model download**. Instantiating `CellposeModel()` fetches the
    `cpsam` weights into the image so runtime tasks never hit the cellpose model
    server. (The v3 `cellpose` component does the same thing via a `docker run`
    line, for the same reason: several concurrent segmentation tasks all fetching
    weights caused intermittent `HTTP 504` from the weight server.)
  - Merges `/src/base/setup_spatialdata_partial.yaml` (spatialdata stack). Unlike
    `cellpose`, it does **not** merge `setup_txsim_partial.yaml` — no txsim needed.
- Nextflow directives: `[ midtime, midcpu, highmem, gpuhighmem ]` — 4 h / 15 CPU /
  100 GB / high-mem GPU. GPU label is what schedules it onto GPU nodes; on CPU the
  SAM model is very slow (the reason v4 was pulled from the CPU test path).

## Arguments (defaults are speed-tuned)

The defaults deliberately trade a little accuracy for speed so the transformer
model finishes in budget:

- `--diameter` (default **30.0**) — expected cell diameter in px. 30 == the model's
  training mean, so it's effectively "no rescaling." Left unset, v4 would run a
  slower size estimate; fixing it at 30 skips that. Bump it for genuinely large
  cells (downsamples → faster + better convergence).
- `--flow_threshold` (default **0.0**) — flow-error QC threshold. Cellpose's own
  default is 0.4. **Setting it to 0 disables the flow-consistency check**, which
  skips a recompute step → faster, at the cost of possibly keeping some ill-shaped
  ROIs. Raise toward 0.4 if you see bad masks.
- `--cellprob_threshold` (default **0.0**) — cell-probability threshold; the
  network emits a per-pixel logit in roughly **−6…+6** and only pixels **above**
  this value seed masks. This is the Cellpose default (0.0). **Lower it** (e.g. −2)
  to recover more / dimmer / low-contrast cells (higher recall); **raise it**
  (e.g. +1…+2) to suppress detections in dim regions (higher precision). Unlike the
  other knobs it's a pure recall↔precision dial, not a speed dial. Exposed here so
  it can be swept per dataset.
- `--niter` (default **10**) — number of flow-dynamics iterations. Cellpose default
  is `None` (≈ proportional to diameter, often ~200). **10 is aggressively low** →
  fast, but pixels for large/elongated cells may not fully converge. Increase if
  cells look under-segmented/split.
- `--min_size` (default **-1**) — minimum pixels per mask. `-1` **disables**
  small-mask removal (Cellpose treats `min_size<0` as "off"); skips a filtering
  pass. Set a positive value (v3 default was 15) to drop specks.
- `--resample` (default **false**) — whether to run dynamics at the original image
  resolution. `false` runs them on the downsampled grid → faster, slightly coarser
  boundaries. `true` gives smoother/more precise masks but is slower.

All six are simply forwarded into `model.eval`. `flow_threshold=0`,
`cellprob_threshold=0` and `min_size=-1` are meaningful values (not "unset") and
*are* passed through — only `None` values get dropped in the eval-params
comprehension (`:54`).

## Optimization / tuning

Two things to understand before tuning:
1. **Every shipped default is biased toward speed**, not accuracy — the component
   exists precisely because the v4 model is slow (see history). So "optimizing for
   quality" here mostly means *walking those defaults back* toward Cellpose's own
   defaults, spending the time budget you can afford on GPU.
2. Grounded in the Cellpose-SAM paper + docs (defaults quoted from the docs). The
   levers below are ranked by expected impact on iST segmentation quality.

**Tier 1 — highest impact on quality**

- **`diameter`** *(exposed; default 30, Cellpose default = auto)*. The biggest
  lever. Cellpose rescales the image so objects land near the ~30 px training mean;
  `30` therefore means "assume objects are already ~30 px / no rescaling." For
  Xenium morphology (~0.2125 µm/px) an ~8 µm nucleus ≈ 35–40 px and a whole cell
  ≈ 60–70 px, so a fixed 30 can be materially wrong. v4 is size-*robust* but not
  size-*invariant* — tune per dataset, or set `None`/`0` to auto-estimate (slower).
- **`cellprob_threshold`** *(now exposed; default 0.0, range −6…+6)*. Pure
  recall↔precision dial. Lower (e.g. −2) → recover more / dimmer / low-contrast
  cells; raise (e.g. +1…+2) → suppress dim-region detections. In iST you usually
  want to capture *all* cells, so sweeping this (≈ −2 … +2) is the most direct way
  to trade missed cells against false ones. No speed cost.
- **`flow_threshold`** *(exposed; default 0.0 here, Cellpose default 0.4)*. `0.0`
  **disables** the flow-consistency QC (fast, keeps ill-shaped masks). Restore
  ~0.4 to filter malformed ROIs; raise above 0.4 to recover more.

**Tier 2 — quality/speed trade-offs (already exposed)**

- **`niter`** *(default 10; Cellpose default None ≈ ∝ diameter, often ~200)*. `10`
  is aggressively low → large/elongated cells may not converge → fragmentation.
  Raise (or set `None`) for non-round cells.
- **`resample`** *(default False)*. `True` runs dynamics at full resolution →
  smoother, more precise boundaries, slower. Turn on when boundary accuracy matters.
- **`min_size`** *(default −1 = off; Cellpose default 15)*. Set positive (≈ 15–50)
  to drop debris/specks → higher precision.

**Tier 3 — not exposed, worth adding for a serious sweep**

- **`normalize`** *(Cellpose default True, 1/99-percentile)*. Accepts a dict
  (percentile bounds, tile-wise normalization, sharpening). For iST images with
  uneven illumination / low contrast, the normalization scheme can matter; expose
  at least the percentiles if results look intensity-sensitive.
- **`batch_size` / `bsize`** *(tile size, 256 for cpsam)*. Throughput knobs — raise
  `batch_size` to trade GPU memory for speed.
- **`augment` / TTA** — test-time augmentation is the accuracy ceiling at ~4–8×
  cost; rarely worth it inside a benchmark.

**Tier 0 — the biggest lever is the input, not a parameter**

The script segments **`image[0]` — a single channel only** (`:56`). Cellpose-SAM
dropped the `channels` arg *because* it's trained to use up to 3 channels
(cytoplasm + nuclear, any order). If a dataset carries a membrane/boundary stain
alongside the nuclear channel, feeding **both** as a 2-channel stack would improve
*whole-cell* boundaries more than any threshold tweak — the current code discards
that signal. This is the highest-leverage change when input images are
multi-channel (but see the "only channel 0" gotcha below — it's not wired up).

**Suggested quality-first sweep order:** `diameter` (per dataset or auto) →
`cellprob_threshold` (−2…+2) → `flow_threshold` → 0.4 → `niter` ↑ / `resample` →
True → `min_size` → ~15–50.

## Risk points / gotchas

- **Speed / GPU.** cpsam is a ViT — on CPU it can blow the `viash test`/CI time
  limit (that's literally why the v3 pin exists on the sibling component). Run this
  one on GPU (`gpuhighmem`). If someone flips the test path to include it on a
  CPU-only runner, expect timeouts.
- **Whole image loaded into RAM** (`:48`) — `.compute().to_numpy()` on `scale0` is
  the full-res plane; big panels are why the label is `highmem`. No tiling/chunking
  is done here beyond what cellpose does internally.
- **Only channel 0 is segmented** (`image[0]`). iST morphology images here are
  effectively single-channel; if a multi-channel image ever arrives, the other
  channels are ignored (v4 could use up to 3 — not wired up).
- **Aggressive defaults can under-segment.** `niter=10`, `flow_threshold=0`,
  `min_size=-1`, `resample=false` are all tuned for throughput. If a dataset gives
  poor masks, the first knobs to relax are `niter` ↑, `flow_threshold` → 0.4,
  `resample` → true.
- **Shared metadata block** (`:67-86`) is duplicated in `cellpose/script.py`; fix
  both together.
- **Don't downgrade to the `Cellpose` class.** It no longer exists in v4; only
  `CellposeModel` is valid. Likewise there is no `channels`/`model_type` in v4.

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: listed as a dependency
  (`methods_segmentation/cellposev4`) and included in the default method string
  `custom_segmentation:cellpose:cellposev4:binning:stardist:watershed`.
- `src/workflows/run_benchmark/main.nf:98`: imported alongside `cellpose`.
- Run scripts: enabled in the GPU/full/MPII/seqera run configs
  (`run_gpu_nebius.sh`, `run_full_nebius.sh`, `run_mpii_nebius.sh`,
  `run_full_seqeracloud.sh`); commented out in `run_test_local.sh` /
  `run_test_nebius.sh` (GPU + time cost too high for the quick test path).

## References

- Docs: https://cellpose.readthedocs.io/en/latest/ (settings, api pages have the
  `model.eval` parameter reference).
- Repo: https://github.com/MouseLand/cellpose
- **Cellpose-SAM (the model this component runs):** Pachitariu M., Rariden M.,
  Stringer C., "Cellpose-SAM: superhuman generalization for cellular
  segmentation", bioRxiv **10.1101/2025.04.28.651001**, v1, 2025-05-01
  (https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1). Not yet a
  peer-reviewed journal version at time of writing.
- **Original Cellpose (the DOI actually in `config.vsh.yaml`):** Stringer et al.
  (2021), Nat. Methods, DOI 10.1038/s41592-020-01018-x — the flow-dynamics CNN,
  i.e. the *framework*, not the v4 model. See the "Citation caveat" up top.
- Adapted from `openproblems-bio/task_spatial_segmentation` (cellpose method).
