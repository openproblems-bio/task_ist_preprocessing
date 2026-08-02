# Cellpose (v3, cyto/nuclei CNN) segmentation — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how it
works, why it goes through the `txsim` wrapper, why it's a *separate* component
from `cellposev4`, what its defaults really are (they are **not** speed-tuned —
see below), and where the tuning levers are.

Links: config `links.documentation`/`links.repository` both point at this task's
own repo (not upstream). Upstream Cellpose: docs
https://cellpose.readthedocs.io/en/latest/ · repo
https://github.com/MouseLand/cellpose.

**Citation caveat (milder than cellposev4's).** The config's `references.doi` is
`10.1038/s41592-020-01018-x` = "Cellpose: a generalist algorithm for cellular
segmentation" (Stringer et al. 2021, Nat. Methods) — the original flow-dynamics
CNN. This component's **default model is `cyto`, which IS that 2021 model**, so the
DOI is appropriate for the shipped default (unlike `cellposev4`, whose DOI pointed
at the wrong paper). Caveat only: the pinned library is Cellpose **3.x**, and the
`--model_type` sweep can select `cyto2` (Cellpose 2.0, Pachitariu & Stringer 2022,
`10.1038/s41592-022-01663-4`) and `cyto3` (Cellpose3, Stringer & Pachitariu 2025,
`10.1038/s41592-025-02595-5`), which come from later papers. If the benchmark
settles on cyto2/cyto3, add those DOIs to the `references` list.

## What this component is

One of the methods at the **segmentation** stage of the iST preprocessing
benchmark (`src/methods_segmentation/`). It implements the standard API
`src/api/comp_method_segmentation.yaml`: `raw_ist.zarr` in → `segmentation.zarr`
out (a 2D cell **label map** + a stub `table`).

It runs **Cellpose v3** (`cellpose<4.0.0`, currently resolves to 3.1.1.x) via the
**txsim wrapper** `txsim.preprocessing.segment_cellpose`, not by calling cellpose
directly. That indirection is the single most important thing to understand about
this component — see the "txsim wrapper" section. The default model is `cyto`
(the original CNN cytoplasm model).

### Why it exists separately from `cellposev4`

There are **two** Cellpose components in this repo, and the split is deliberate.
(This table is the mirror of the one in `cellposev4/NOTES.md`.)

| | `cellpose/` (this one) | `cellposev4/` |
|---|---|---|
| cellpose version | `cellpose<4.0.0` (pinned) | `cellpose>=4.0.0` |
| model | `cyto` (CNN) via `models.Cellpose` | `cpsam` (ViT/SAM) via `CellposeModel` |
| call path | `txsim.preprocessing.segment_cellpose(...)` wrapper | direct `model.eval(...)` |
| extra deps | needs `txsim` | no `txsim` (leaner image) |
| arg surface | ~18 args (channels, model_type, invert, do_3D, …) | 6 args (diameter, flow_threshold, cellprob_threshold, niter, min_size, resample) |
| defaults | == Cellpose library defaults (quality) | deliberately **speed-tuned** |

History (git `-- src/methods_segmentation/cellpose/`):
- `5487edcfa` "cellpose added" (2024-11) — original component; ~17 args, `base_python`,
  txsim wrapper, `models.Cellpose` (v3 API).
- `35061648c` "Support cellpose v4" (2025-06) — reshuffled args to be v4-compatible
  (moved `resample`, commented out `interp`/`net_avg`, note "diameter should be None
  with v4"). This is why some args are commented out in the config today.
- `8015045ab` "Use cellpose v3 instead of v4 to stay in test time limits" — pinned
  `cellpose<4.0.0`, because the v4 SAM/ViT model was **too slow** to finish inside
  the CI/`viash test` budget. This component stayed v3; `cellposev4` (#166) was added
  as the heavier GPU sibling.
- `7e850cc84` "Add gpu for cellpose" + `702ba3fd2` — added the `gpu` label and moved
  to `base_pytorch_nvidia:1`. **But see the GPU gotcha below** — the txsim v3 path
  does not actually use the GPU.
- `31d57ec1c` "Add method selection arguments" — merged `setup_spatialdata_partial`.

## How Cellpose v3 works (background)

Cellpose predicts, per pixel, a **flow vector** toward the center of the cell that
pixel belongs to, plus a **cell probability**. At inference each pixel "follows the
flow" for `niter` iterations to a fixed point; pixels converging to the same point
become one **mask**. This flow-dynamics core is shared across Cellpose 1–4.

v3 exposes two model classes (both still present in `cellpose<4`):
- **`models.Cellpose`** — bundles a `SizeModel` (auto-diameter estimator) with a
  `CellposeModel`. Constructor default `gpu=False`, `model_type="cyto3"`. This is
  what the txsim wrapper uses.
- **`models.CellposeModel`** — the mask model alone (what `cellposev4` uses).

Builtin models (`MODEL_NAMES`, cellpose 3.1.1.1): `cyto3`, `nuclei`, `cyto2`,
`cyto` (+ specialist `*_cp3`, tissuenet, livecell, bacterial, transformer models).
`cyto`/`nuclei`/`cyto2`/`cyto3` each have a paired size model, so auto-diameter
works for them. `diam_mean` is 30 px for cyto* models, 17 px for `nuclei`.

## txsim `segment_cellpose` wrapper — the key indirection

`script.py` calls `tx.preprocessing.segment_cellpose(image[0], hyperparameters)`.
The wrapper lives in `theislab/txsim@dev`
(`txsim/preprocessing/_segmentation.py`, function `segment_cellpose`). For
`cellpose<4` it does, in order:

1. Reads `hyperparams["model_type"]` (defaults to `'nuclei'` if absent — but this
   component always passes `cyto`).
2. `model = models.Cellpose(model_type=model_type)` — **note `gpu=` is NOT passed**,
   so it defaults to `gpu=False`. See the GPU gotcha.
3. `del hyperparams["model_type"]` — the model_type key is consumed here and does
   **not** reach `eval`.
4. `res, _, _, _ = model.eval(img, channels=[0, 0], **hyperparams)` —
   **`channels=[0,0]` is hardcoded** (segment `img` as a single grayscale plane, no
   separate nuclear channel). Every remaining config arg is forwarded as a kwarg.

So the contract is: **all config args except `model_type` are forwarded verbatim
into `Cellpose.eval`** (which passes the non-explicit ones through `**kwargs` into
`CellposeModel.eval`). Adding a new config arg whose name matches an `eval` kwarg
is therefore enough to wire it — no script change needed (that is exactly how the
newly added `--niter` works).

## script.py — step by step

1. **Build hyperparameters** (`:34-38`) — `hyperparameters = par.copy()`, then
   `{k:(v if v != "None" else None)}` converts the string `"None"` (how the
   string-typed args `channel_axis`, `z_axis`, `rescale`, `anisotropy` carry "unset")
   back to Python `None`, then `del`s `input`/`output`. Everything else — all the
   segmentation knobs — stays in the dict.
2. **Read input** (`:40-43`) — `sd.read_zarr(par["input"])`;
   `sdata['image']['scale0'].image.compute().to_numpy()` pulls the full-res image
   into RAM; the `scale0` transform is copied to reattach to the output.
3. **Segment** (`:44`) — `tx.preprocessing.segment_cellpose(image[0], hyperparameters)`.
   `image[0]` = first channel of a `(c, y, x)` array → a single 2D plane. Returns
   the label array.
4. **Downcast + parse** (`:45-48`) — `convert_to_lower_dtype` shrinks the label array
   to the smallest uint that holds `max label`; wrap as `xarray.DataArray(dims=('y','x'))`,
   `Labels2DModel.parse` with the copied transform, store as
   `sd_output.labels['segmentation']`.
5. **Carry the metadata table** (`:50-69`) — copies a stub `table` from the input's
   `tables["metadata"]`. Downstream needs a `cell_id` column; resolved in priority
   order: explicit `obs["cell_id"]` → the table's `instance_key` column
   (`uns["spatialdata_attrs"]["instance_key"]`) → the obs index. Optional `region`
   column carried if present. **This block is duplicated verbatim in
   `cellposev4/script.py`** — the instance_key fallback was added for exports (e.g.
   the Xenium WTA preview behind the Atera dataset) that lack an explicit `cell_id`.
   Keep the two copies in sync.
6. **Write** (`:71-74`) — `rmtree` any existing output, then `sd_output.write`.

Note: the output `table` describes the *reference* cells from metadata, not the
newly predicted masks — segmentation output here is the **label image**; the
counting/aggregation stage downstream turns masks into a cell × gene table.

## config.vsh.yaml — the Docker setup

- Base image `openproblems/base_pytorch_nvidia:1` (CUDA/torch), GPU-capable — but
  the v3 code path does not use the GPU (gotcha below).
- `setup`:
  - `pypi: cellpose<4.0.0` — the version pin that makes this "v3".
  - `docker.run`: a **build-time model download** — instantiates
    `models.Cellpose(gpu=False, model_type=m)` for `m in ['cyto','nuclei','cyto2','cyto3']`,
    which caches each mask model **and its size model** into the image. This was
    extended (this tuning session) from caching only `cyto` to caching all four,
    because the `--model_type` sweep selects among them and several concurrent tasks
    fetching weights caused intermittent `HTTP 504` from the weight server.
  - Merges `/src/base/setup_txsim_partial.yaml` (txsim + squidpy + rasterio) **and**
    `/src/base/setup_spatialdata_partial.yaml`. Unlike `cellposev4`, this one needs
    txsim (for the wrapper).
- Nextflow directives: `[ midtime, midcpu, highmem, gpuhighmem ]` — 4 h / 15 CPU /
  100 GB / high-mem GPU. The `gpuhighmem` label schedules onto a GPU node, but see
  the GPU gotcha — the GPU sits idle for the v3 path.

## Arguments

The forwarded `Cellpose.eval` / `CellposeModel.eval` args, their config defaults,
and how they compare to Cellpose 3.1.1.1's own defaults. **Almost every default
here already equals Cellpose's library default** — the only deviation is
`model_type` (`cyto` vs library `cyto3`). So, unlike `cellposev4`, this component
is **not** shipped with speed-tuned defaults; it is essentially "stock Cellpose."

| arg | config default | Cellpose default | notes |
|---|---|---|---|
| `--model_type` | `cyto` | `cyto3` | selects the model; consumed by the wrapper, not passed to eval. `cyto`/`nuclei`/`cyto2`/`cyto3`. **deviation from library default** |
| `--diameter` | `30.0` | `30.` (class) | fixed size assumption. `0`/`None` → auto-estimate via SizeModel (slower). biggest lever |
| `--flow_threshold` | `0.4` | `0.4` | flow-error QC; higher = keep more (more permissive), lower = stricter shape filter |
| `--cellprob_threshold` | `0.0` | `0.0` | per-pixel logit gate (~−6…+6); lower = more/dimmer cells (recall), higher = fewer (precision) |
| `--min_size` | `15` | `15` | drop ROIs below N px; `0` keeps specks |
| `--resample` | `True` | `True` | run dynamics at full resolution (smoother, slower). already the quality setting |
| `--normalize` | `True` | `True` | 1/99-percentile normalization per channel |
| `--niter` | `0` (auto) | `None` (auto) | **newly exposed.** dynamics iterations; 0/None → ∝ diameter (~200 at resample). lower for speed |
| `--augment` | `False` | `False` | test-time augmentation (TTA); quality at ~4–8× cost |
| `--batch_size` | `8` | `8` | tiles per forward pass |
| `--tile_overlap` | `0.1` | `0.1` | tile overlap fraction |
| `--invert` | `False` | `False` | invert intensities before running |
| `--rescale` | `None` | `None` | manual resize factor (only used if diameter is None) |
| `--do_3D` | `False` | `False` | 3D segmentation (input here is 2D) |
| `--anisotropy` | `None` | `None` | 3D only |
| `--stitch_threshold` | `0.0` | `0.0` | 3D stitching of 2D masks |
| `--channel_axis` | `None` | `None` | auto-detected; input is a single 2D plane |
| `--z_axis` | `None` | `None` | 3D only |

Commented-out in the config: `--net_avg` (removed in cellpose 2.2+), `--tile`
(deprecated), `--interp` (valid in v3 — default `True` — but commented out during
the abandoned v4-compat pass; low sweep value since flipping it off only degrades
2D dynamics).

## Optimization / tuning

Two things to understand before tuning (both differ from the `cellposev4` story):

1. **Defaults are already Cellpose's own defaults, i.e. quality-oriented** (flow QC
   on at 0.4, min-mask filtering on at 15, resample on, normalize on). The only
   non-stock default is `model_type=cyto`. So "optimize for quality" here is **not**
   "walk speed-tuned defaults back to library defaults" (that was the v4 job) — it's
   *exploring the model choice, the object scale, and the recall/precision dials*.
2. Grounded in the cellpose 3.1.1.1 source (`models.py`, `dynamics.py`) and docs,
   not memory. Ranges tied to Xenium morphology (~0.2125 µm/px).

**Tier 0 — the input, not a parameter.** The wrapper hardcodes `channels=[0,0]` and
the script feeds only `image[0]` — a **single grayscale plane**. iST morphology is
effectively single-channel here, but if a dataset carried a membrane/boundary stain
alongside the nuclear one, feeding both (via `channels=[cyto_ch, nuc_ch]`) would
improve whole-cell boundaries more than any threshold tweak. Not wired up (would
require changing the hardcoded `channels` in txsim, which is out of this repo).

**Tier 1 — highest impact on quality**

- **`model_type`** *(default `cyto`; library default `cyto3`)*. The most
  distinctive v3 lever (v4 has no model choice). On a single grayscale morphology
  channel the model matters a lot: `nuclei` for a nuclear-dominant stain, `cyto2`
  (Cellpose 2.0) and `cyto3` (Cellpose3 generalist super-model) as improved
  generalists over the 2021 `cyto`. Sweep `nuclei`/`cyto2`/`cyto3`.
- **`diameter`** *(default 30.0; `0`/`None` = auto)*. Cellpose rescales so objects
  land near the model's `diam_mean` (30 px cyto, 17 px nuclei); a fixed 30 assumes
  ~30 px objects. Xenium (~0.2125 µm/px): an ~8 µm nucleus ≈ 38 px, a whole cell
  ≈ 60–70 px, so 30 can be materially wrong. Sweep `0` (auto), ~40 (nucleus scale),
  ~60 (whole-cell scale).
- **`cellprob_threshold`** *(default 0.0, range −6…+6)*. Pure recall↔precision dial,
  no speed cost. Lower (−1…−2) recovers more/dimmer cells; raise (+1…+2) suppresses
  dim detections. In iST you usually want all cells → this is the most direct dial.
- **`flow_threshold`** *(default 0.4)*. Already at the Cellpose default (QC on).
  Raising it (0.6/0.8) keeps cells with higher flow error (more, possibly
  ill-shaped); lowering it (0.2) is a stricter shape filter.

**Tier 2 — quality/speed trade-offs**

- **`niter`** *(newly exposed; default 0 = auto ≈ 200 at resample)*. Lower (50) is
  faster — meaningful because the v3 path is CPU-bound (see gotcha); raise for
  large/elongated cells that under-converge.
- **`resample`** *(default True = quality setting)*. Only non-default is `False`
  (dynamics on the downsampled grid → faster, coarser boundaries). Include to
  characterize the speed cost of the current default.
- **`min_size`** *(default 15)*. `0` keeps small specks (higher recall on tiny
  cells); larger (50+) drops debris (higher precision).
- **`augment`** *(default False)*. TTA — the accuracy ceiling at ~4–8× cost; rarely
  worth it in a benchmark, include one `True` to bound the gain.
- **`normalize`** *(default True)*. Flip to `False` to see intensity sensitivity;
  usually worse for uneven-illumination iST images.

**Tier 3 — not exposed, worth adding for a serious sweep**

- **`normalize` as a dict** *(percentiles, `tile_norm`, `sharpen`)*. cellpose v3
  accepts a normalization dict (see `models.normalize_default`). For iST images with
  uneven illumination, tile-wise normalization / sharpening can matter; the boolean
  `--normalize` can't express it. Would need a JSON/dict-typed arg.
- **`bsize`** *(224)* / **`max_size_fraction`** *(0.4)* — throughput / large-mask
  culling; low value for a first sweep.

**Suggested quality-first order:** `model_type` (nuclei/cyto2/cyto3) → `diameter`
(auto or data scale) → `cellprob_threshold` (−2…+2) → `flow_threshold` →
`min_size`. The current sweep (`scripts/run_benchmark/cellpose_params.yaml`)
covers all of these; **20 variants** total (1 default + 19).

## Risk points / gotchas

- **GPU is scheduled but NOT used (v3 path).** The txsim wrapper calls
  `models.Cellpose(model_type=...)` without `gpu=`, so it defaults to `gpu=False` and
  runs on **CPU** — even though the config uses `base_pytorch_nvidia` and the
  `gpuhighmem`/`gpu` label puts it on a GPU node. As of `txsim@dev` at time of
  writing, this wastes the GPU allocation and makes runs slow. Fixing it properly
  means either patching txsim to pass `gpu=core.use_gpu()` (as the v4 branch of the
  same wrapper does) or dropping the GPU label. Left for the user — do not silently
  flip infra labels. This is also why the `--niter` speed knob matters here.
- **`--model_type` beyond `cyto` triggers weight downloads unless pre-cached.** The
  image now pre-caches `cyto`/`nuclei`/`cyto2`/`cyto3` (config `docker.run`), so the
  sweep is safe; if you add another model to the sweep, add it to that list too or
  concurrent tasks may hit `HTTP 504` from the weight server.
- **Whole image loaded into RAM** (`:42`) — `.compute().to_numpy()` on `scale0` is
  the full-res plane; big panels are why the label is `highmem`.
- **Only channel 0 is segmented, as grayscale** (`image[0]` + hardcoded
  `channels=[0,0]`). Any additional stain channels are ignored.
- **`model_type` is consumed by the wrapper, not passed to eval.** If you rename or
  drop it, the wrapper falls back to `'nuclei'` — a silent behaviour change.
- **String "None" convention.** `channel_axis`/`z_axis`/`rescale`/`anisotropy` are
  typed `string` with default `"None"`; the script converts the literal `"None"` →
  Python `None`. A new numeric knob that needs an "auto/unset" sentinel should reuse
  a real numeric sentinel the tool understands (as `--niter 0` does), not the string
  hack, or it will be forwarded as a string and break eval.
- **Shared metadata block** (`:50-69`) is duplicated in `cellposev4/script.py`; fix
  both together.
- **Don't confuse the classes.** This component uses `models.Cellpose` (with
  SizeModel); `cellposev4` uses `models.CellposeModel`. `Cellpose` still exists in
  v3 but was removed in v4.

## Wiring

- `src/workflows/run_benchmark/config.vsh.yaml`: listed as a dependency
  (`methods_segmentation/cellpose`, L158) and included in the default segmentation
  method string `custom_segmentation:cellpose:cellposev4:binning:stardist:watershed`
  (L65).
- `src/workflows/run_benchmark/main.nf:97`: imported alongside `cellposev4`.
- Run scripts: enabled in `run_full_local.sh`, `run_full_nebius.sh`,
  `run_full_seqeracloud.sh`, `run_gpu_nebius.sh`, `run_mpii_nebius.sh`,
  `run_test_seqeracloud.sh`; commented out in `run_test_local.sh` /
  `run_test_nebius.sh`. Dedicated sweep runners added this session:
  `scripts/run_benchmark/run_test_cellpose_local.sh` and `_nebius.sh`, driven by the
  committed `scripts/run_benchmark/cellpose_params.yaml`.

## References

- Docs: https://cellpose.readthedocs.io/en/latest/ (settings + api pages have the
  `model.eval` parameter reference).
- Repo: https://github.com/MouseLand/cellpose (defaults verified against tag
  `v3.1.1.1`, `cellpose/models.py` + `cellpose/dynamics.py`).
- txsim wrapper: `theislab/txsim@dev`, `txsim/preprocessing/_segmentation.py`,
  `segment_cellpose`.
- **Original Cellpose (`cyto`, the config DOI):** Stringer, Wang, Michaelos,
  Pachitariu (2021), Nat. Methods, **10.1038/s41592-020-01018-x**.
- **Cellpose 2.0 (`cyto2`):** Pachitariu & Stringer (2022), Nat. Methods,
  **10.1038/s41592-022-01663-4**.
- **Cellpose3 (`cyto3`, image restoration):** Stringer & Pachitariu (2025), Nat.
  Methods, **10.1038/s41592-025-02595-5**.
