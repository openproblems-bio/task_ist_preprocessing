# Segger transcript-assignment — implementation notes

Reference for anyone (human or Claude) working on this component. Captures how it
works, why the Docker setup looks the way it does, and where it can break.

## What this component is

One of 7 methods at the **transcript-assignment** stage of the iST preprocessing
benchmark (`src/methods_transcript_assignment/`). It implements the standard API
`src/api/comp_method_transcript_assignment.yaml`: raw iST + segmentation in →
`transcript_assignments.zarr` out (see `src/api/file_transcript_assignments.yaml`).

Unlike the sibling methods, segger is a **GPU-only, GNN-based** approach. This
component is an **adapter** around segger's own end-to-end CLI (`segger segment`):

- The *segmentation* stage produces a label map (cell boundaries).
- Segger's job here is to **re-assign transcripts to cells**, using that
  segmentation only as a **boundary prior**, then letting its GNN refine.
- Flow: standard SpatialData `.zarr` in → synthesize a **Xenium-layout directory**
  → run `segger segment` as a subprocess → map results back onto the standard
  `transcript_assignments.zarr`.

Links: docs https://elihei2.github.io/segger_dev/ · repo https://github.com/dpeerlab/segger
· DOI 10.1101/2025.03.14.643160

## script.py — step by step

0. **Defensive `CUDA_VISIBLE_DEVICES` guard** (`:44-51`) — if CVD is set-but-**empty**,
   delete it (empty masks all GPUs from the CUDA driver API while torch's NVML gate still
   passes). **This is NOT the observed GPU bug** (the real Nextflow task runs with CVD
   *unset*, verified via the failed task's env) — it's cheap insurance against launchers that
   do pass an empty value. The real GPU fix is the `NUMBA_CUDA_USE_NVIDIA_BINDING` env var in
   step 4; see the gotcha below.

1. **GPU gate** (`:56-62`) — raises immediately if `torch.cuda.is_available()` is
   false. Segger needs the GPU end-to-end (cudf, cuspatial, torch kernels); fail
   fast instead of crashing deep in the subprocess. `viash test` therefore cannot
   run on a CPU-only host. **Caveat:** this gate is NVML-based, so it does NOT catch a GPU
   the CUDA *driver* API can't use (empty CVD, or the numba-binding issue in step 4).

2. **Boundaries → `{cell,nucleus}_boundaries.parquet`** (`:63-93`) —
   `sd.to_polygons(segmentation)`, transformed into the transcripts' coordinate
   frame via a composed `Sequence([segm→global, global→transcripts.inverse()])`,
   written in Xenium's **per-vertex** schema (`cell_id` **string**, `vertex_x/y` f64).
   `cell_id` MUST be a string (cast via int64 → `"1"`/`"2"`) to match the transcripts'
   string `cell_id` (step 3): segger joins transcript `cell_id` → boundary `cell_id` to
   attach each cell to its polygon (`setup_anndata`), and a str-vs-int mismatch makes every
   row miss → all-`nan` `entity_index` → `Index has duplicate keys: [nan]`.
   No cell/nucleus split exists, so the **same polygon set is written to both**
   files; `--prediction_mode` selects which drives the graph.

3. **Transcripts → `transcripts.parquet`** (`:95-153`) — writes a Xenium-schema
   transcript table: `transcript_id`, `x/y[/z]_location`, `feature_name`,
   `cell_id` (prior; background → `"UNASSIGNED"`), plus dummy `qv=40.0` and
   `overlaps_nucleus`. The **prior `cell_id`** per transcript comes from looking up
   the label image at each transcript's integer pixel coords (truncated to int64,
   then **clamped** to image bounds — a few can round a pixel past the raster edge, same
   clamp as `basic_transcript_assignment` / `baysor` / `proseg`). `tx_pd` (clean RangeIndex)
   is the **canonical transcript frame** over ALL transcripts, and **ALL of them are written
   to the parquet**, so the parquet row order == `tx_pd` and segger's reported `row_index`
   indexes it directly (step 5). Guards against a total coordinate-frame mismatch by raising
   if *every* transcript clamps out of bounds (`n_oob == n_tx`).
   **History:** this step formerly EXCLUDED out-of-bounds transcripts from the parquet and
   remapped `row_index` through a `seg_orig_idx` array — a workaround for when the
   process_dataset crop left ~38% of transcripts outside the cropped labels. That crop bug is
   fixed (`crop_points_by_global_xy` crops transcripts to the same global box as the labels)
   and a full run on re-processed Xenium + MERFISH confirmed `n_oob == 0`, so the exclusion
   was a no-op and was removed in favour of the plain clamp. This is INDEPENDENT of the
   empty-`bd`-batch crash (see that gotcha) — that one is intra-field, not OOB-driven.

4. **Run segger** — `segger segment -i -o --n-epochs --prediction-graph-buffer-ratio
   --prediction-mode --node-representation-dim`, launched via a `run_segger(node_dim)`
   helper (`Popen`, stderr→stdout, **streamed live AND captured**). **Sets
   `NUMBA_CUDA_USE_NVIDIA_BINDING=1` — load-bearing (see gotcha): without it numba
   enumerates 0 GPUs inside cudf's `.values`/`as_cuda_array` path and crashes in
   `phenograph_rapids`.** Also sets
   `RAPIDS_NO_INITIALIZE / CUDF_NO_INITIALIZE / RMM_NO_INITIALIZE = 1` so
   `import cudf` (eager during tiling) doesn't abort on driver validation.
   **PCA-ceiling auto-retry:** segger runs `PCA(n_components=node_dim)` on a
   **genes × genes** correlation matrix in `setup_anndata`. On small panels/crops fewer
   genes survive segger's internal count filters (`cells_min_counts=10`, gene count filter)
   than `node_dim`, so the PCA raises `n_components=128 must be between 0 and
   min(n_samples, n_features)=<G>`. We can't cheaply predict `<G>` from here (depends on
   segger's own filtering + order), but segger prints the exact ceiling — so on a non-zero
   exit we regex `min(n_samples, n_features)=(\d+)` out of the captured output and **retry
   once** at that dim. On full panels (≥`node_dim` genes survive) this never fires. Then
   asserts `segger_segmentation.parquet` exists.

5. **Map back** (`:181-235`) — read segger output, filter to confident assignments:
   **`segger_similarity >= similarity_threshold`** (segger 0.1.0 dropped the old
   `keep` column; this is the exact filter segger itself applies when writing its
   AnnData) AND drop null / `-1` / `UNASSIGNED` / `NONE` / `NAN` / `""` cell ids.
   Use segger's `row_index` column to scatter `pd.factorize(segger_cell_id,
   sort=True)` codes **`+1`** back onto `tx_pd` (background stays `0`). Build output
   `SpatialData`: `points["transcripts"]` (native coords + `cell_id`) +
   `tables["table"]` AnnData whose `obs` lists the unique assigned cell IDs with `0`
   excluded. Temp working dir is removed at end.

### Duplicate-index workaround (shared with basic_transcript_assignment)

`:102-111`. Multi-partition parquet each start a 0-based index → duplicate indices
in the combined dask frame → `sd.transform` fails ("cannot reindex on an axis with
duplicate labels"). Fix: `.compute().reset_index(drop=True)` then rebuild as a
single-partition dask frame with a clean RangeIndex before transforming. Do **not**
transform `sdata[transcripts_key]` directly.

## Arguments

| arg | default | notes |
|-----|---------|-------|
| `--transcripts_key` | `transcripts` | key in `sdata.points` |
| `--coordinate_system` | `global` | must exist on both transcripts and segmentation |
| `--n_epochs` | `20` | GNN training epochs |
| `--prediction_graph_buffer_ratio` | `0.05` | fraction of each polygon's equiv. radius to buffer/expand at prediction → segger's `--prediction-graph-buffer-ratio` (renamed from `--prediction-expansion-ratio`; default was 0.5) |
| `--prediction_mode` | `cell` | one of `nucleus` / `cell` / `uniform` (matches segger's own 0.1.0 default) |
| `--node_representation_dim` | `128` | segger's `--node-representation-dim` (== `in_channels`/`cells_embedding_size`); PCA n_components for the gene/cell embeddings. Effectively an **upper bound** — script.py auto-retries at segger's reported ceiling on small panels. Leave 128 for full panels. |

## config.vsh.yaml — the Docker setup is load-bearing

segger is fundamentally a **RAPIDS application** (needs cudf + cuspatial + cugraph +
cuml). The image is therefore based on the **official RAPIDS image**, not NGC PyTorch.
**Order matters — don't reorder the setup steps, and don't add `apt` steps** (see below).

### Why RAPIDS base, not NGC PyTorch (history)

The image formerly bolted the RAPIDS stack onto `nvcr.io/nvidia/pytorch:26.06`. That got
far but hit a **wall**: segger's `cugraph → raft_dask` import loads UCX, which collided
with the NGC image's bundled HPC-X UCX (`/opt/hpcx/ucx`) and aborted with heap corruption
(`free(): invalid pointer`, **SIGABRT / exit 134**). No pip pin fixes a native-library
conflict. Rebasing on `rapidsai/base` (which ships a consistent, working UCX) resolves it
— validated: the exact `cugraph.dask.comms` import that crashed now exits 0.

### The recipe

- **Base `rapidsai/base:25.04-cuda12.8-py3.12`** — brings cudf/cugraph/cuml/cuspatial
  25.04 + working UCX/dask, plus **numpy 2.0.2 / pandas 2.2.3 / python 3.12** for free.
- **⚠️ Runs as non-root `rapids` (uid 1001)** — so **`type: apt` does NOT work** in a
  viash build (needs root; the build inherits the base's USER). Everything installs via
  **conda / pip only** (both write the rapids-owned `/opt/conda`). Consequences:
  - **git** (needed for the segger + openproblems github pip installs) →
    `conda install -c conda-forge git`, not apt.
  - **cv2 / libGL**: segger pulls `opencv-python`, whose `cv2` needs `libGL.so.1` (a
    system lib we can't apt-install). Fix: swap to **`opencv-python-headless`** (same
    `cv2`, no libGL) as the last step.
  - **⚠️ RUNTIME: can't write Nextflow's scratch dir.** As uid 1001 the container can't
    write the root-owned task work dir → `Permission denied` on `.command.log` /
    `.command.begin` / `.exitcode` (fails before `script.py` runs). Since viash can't set
    the image `USER`, this is fixed **outside** this config: `labels_nebius.config` sets
    `[runAsUser: 0]` on the pod directive of the GPU label segger runs on. **segger currently
    uses the `gpuh100` label** (moved off `gpu` in `58912e46`/`e9505f13a` for H100 memory), so
    `runAsUser: 0` is carried on **both `gpu` and `gpuh100`** (+ `gpuhighmem` as insurance).
    Moving segger to any `gpu*` label without this directive reintroduces the `Permission
    denied` on `.command.log` / `.command.begin` / `.exitcode`. Carry `runAsUser: 0` with it.
  - **⚠️ RUNTIME: `/dev/shm` too small for the training DataLoader.** segger's Lightning
    training uses a **torch_geometric multi-worker DataLoader** that hands each collated graph
    `Batch` to the main process through **`/dev/shm`** (torch's file-descriptor sharing
    strategy). Kubernetes defaults a container's `/dev/shm` to **64 MB**, which the first batch
    overflows → `RuntimeError: unable to write to file </torch_...>: No space left on device
    (28)` / `Unexpected bus error … insufficient shared memory (shm)`, killing the DataLoader
    workers in the pre-train **sanity check** (right after the model is built — GPU already in
    use, so it looks like training "just died"). Fixed **outside this config**, same place as
    `runAsUser`: `labels_nebius.config` mounts a **Memory-backed `emptyDir` at `/dev/shm`**
    (`[emptyDir: [medium: 'Memory', sizeLimit: '16Gi'], mountPath: '/dev/shm']`) on the
    `gpuh100` pod directive. **Note the `*sharedmem` labels' `--shm-size` `containerOptions`
    do NOT help** — `containerOptions` is a docker/singularity concept the **k8s executor
    ignores**; on k8s `/dev/shm` must be enlarged via a pod volume. Carry this with
    `runAsUser: 0` if segger changes labels. Fallback if 16Gi ever proves too small on a huge
    dataset: bump `sizeLimit`, or make segger's DataLoader use `num_workers=0` / the
    `file_system` sharing strategy (needs an image change; the shm bump does not).
- **torch 2.6.0 + torchvision 0.21.0 (cu124)** from the PyTorch index — *standard*
  wheels (numpy-2 compatible) that coexist with the RAPIDS cuda-12 libs; torchvision must
  match torch (segger's `data_module` imports it). Big simplification vs NGC: **`torch_scatter`
  has a PREBUILT wheel** (`data.pyg.org/whl/torch-2.6.0+cu124.html`) — **no CUDA compile**.
  Plus `torch_geometric` + `lightning` + `typer` (pure Python — segger imports all three but
  declares none). These GNN-side deps only surface **past** the GPU-gated `rmm.reinitialize`
  in segger's `__init__`; validate them on a CPU pod by mocking `rmm.reinitialize` so the
  full `segger.data` / `segger.cli.segment` import chain runs.
- **spatialdata pinned `<0.8` (→ 0.7.1)** — the dask tension: RAPIDS 25.04 hard-pins
  **`dask==2025.2.0`**, but spatialdata **0.7.2+ require dask ≥2026.x** (only 0.7.0/0.7.1
  accept 2025.2.0). So pin spatialdata `<0.8` AND `dask==2025.2.0`/`distributed==2025.2.0`
  so nothing drags dask up (which breaks dask-cudf/cugraph). `sd.transform` validated on 2025.2.0.
- **segger** pinned to commit `dpeerlab/segger@0233cf62` (== 0.1.0). Formerly unpinned
  trunk, which broke the component when 0.1.0 renamed CLI flags + changed the output
  schema. `script.py` targets THIS commit — re-check `segger segment --help` and
  `data/writer.py` before bumping.
- **pandas / anndata / scanpy pinned as the VERY LAST step** — the same cudf-rooted
  three-way conflict: cudf 25.4 needs **`pandas<2.2.4`**, but anndata ≥0.13 / scanpy ≥1.12
  (pulled by segger) need **`pandas>=2.3`** (anndata 0.13's zarr reader calls
  `pd.StringDtype(na_value=...)`, a pandas-2.3 API → `StringDtype ... unexpected keyword
  'na_value'` reading any zarr table). Hold anndata at **0.12.x** (`pandas>=2.1`) and
  scanpy at **1.11.x** (`pandas>=1.5`). Do it LAST so it wins over segger's unpinned deps.

**Validated coexisting set** (imports on the RAPIDS base, exit 0): pandas 2.2.3, dask
2025.2.0, spatialdata 0.7.1, anndata 0.12.19, scanpy 1.11.5, torch 2.6.0+cu124, torchvision
0.21.0, torch_scatter 2.1.2, typer, cudf/cugraph/cuml/cuspatial 25.04, segger 0.1.0. The full
`segger.data`/`segger.cli.segment` import chain passes with `rmm.reinitialize` mocked. Not yet
run end-to-end on a GPU (rmm/cuInit and the actual segmentation need a GPU node). **Watch:**
numpy is unpinned — the RAPIDS base ships 2.0.2 but the CI-built image drifted to 2.4.6; it
imports fine but cudf 25.4's tested numpy range is lower, so pin numpy if runtime issues appear.

**Cross-version caveat:** the rest of the benchmark uses spatialdata 0.8.0 (dask 2026.x);
segger reads/writes zarr with 0.7.1. SpatialData format versioning should bridge this, but
it's the main thing to watch on the first real run.

Engines: `docker` + `native` (native lets the script run directly on a GPU host).
Nextflow runner labels: `[hightime, midcpu, highmem, gpuh100]`.

## Wiring

- Registered in `src/workflows/run_benchmark/config.vsh.yaml` (dep list `:163`,
  default methods string `:71`) and in the assignment fan-out `main.nf:142`.
- Test/run script: `scripts/run_benchmark/run_test_segger_nebius.sh` — default
  pipeline + segger at the assignment stage, on **Nebius GPU** nodes via Seqera
  (`tw launch`, `--config src/base/labels_nebius.config`), against S3 test resources.
- If a config/script edit "doesn't take effect" on the cluster, the image is likely
  stale — use the `check-component` skill (committed to origin/main → regenerated in
  build/main → rebuilt on ghcr `build_main` tag).

## Risk points / gotchas

- **Never start a continuation line in `script.py` with `|` — viash's Nextflow codegen
  deletes it (Groovy `stripMargin`).** viash embeds `script.py` into the Nextflow module
  (`target/nextflow/.../segger/main.nf`) as a Groovy string and runs it through
  `.stripMargin()`, whose **default margin delimiter is `|`**: any line matching
  `^\s*\|` has the leading whitespace **and the `|`** stripped. A wrapped boolean/bitwise
  expression like
  `n_oob = int(np.count_nonzero(\n    (y_coords < 0) | (y_coords >= H)\n    | (x_coords < 0) | (x_coords >= W)\n))`
  becomes `... (y_coords >= H)\n (x_coords < 0) ...` in the generated module — the two
  parenthesised groups then **juxtapose into a call**, so at run time Python raises
  `TypeError: 'numpy.ndarray' object is not callable` (the caret points at the *first*
  group, the "callable"). **This corrupts ONLY the `=nextflow=>` target**, so `viash test` /
  `viash run` (the `=executable=>` target, no `stripMargin`) pass clean and the bug is
  invisible until a Nextflow/Tower benchmark run. It also survives across builds: the
  `build/main` executable can be correct while the `build/main` *nextflow module* is
  corrupt. **Fix = put the operator at the END of the line** (`... (y_coords >= H) |` /
  next line `(x_coords < 0) ...`) so no line begins with `|`; or keep it on one line.
  Grep guard before committing: `grep -nE '^[[:space:]]*\|' script.py` must be empty.
  (Diagnosed 2026-08-06: reproduced byte-for-byte, incl. the exact caret placement.)
- **numba must use the NVIDIA (cuda-python) driver binding, or cudf/numba see 0 GPUs.**
  The load-bearing GPU fix: `run_segger` sets **`NUMBA_CUDA_USE_NVIDIA_BINDING=1`** for the
  subprocess. RAPIDS (cudf/cuml/cugraph) builds its CUDA context through the cuda-python
  binding; numba defaults to its *own* ctypes driver, which then enumerates **zero** devices
  inside cudf's `.values` → `data_array_view` → `numba.cuda.as_cuda_array` path →
  `IndexError: list index out of range` at `numba_cuda …/devices.py self.gpus[devnum]`.
  The crash surfaces in segger's `setup_anndata` → `phenograph_rapids` (`neighbors.py`,
  `result.sort_values('vertex')['partition'].values.get()`) — the first heavy cudf op — even
  though **torch, `cudf.sum()` (libcudf path) and cupy all see the GPU fine**, and even
  though bare `from numba import cuda; cuda.gpus` enumerates 1. Diagnosed live on a GPU pod by
  bisection: not root/`runAsUser`, not the node, not CVD (unset in the real task), not the
  no-init vars, not the numba/numpy versions — only `NUMBA_CUDA_USE_NVIDIA_BINDING` toggles
  it (`unset`/`0` → crash; `1` → cudf `.values.get()` + cuml kNN + cugraph louvain all pass).
  A benign `TypeError: cuDriverGetVersion() takes no arguments (1 given)` / "Not patching
  Numba" warning prints with the binding on — ignore it. This is what the earlier (mistaken)
  "empty CVD" note was actually about; step 0's CVD guard is now just defensive.
- **numba/numpy drift out of the RAPIDS-supported range (pinned in config).** scanpy's
  umap-learn/pynndescent (no numba upper bound) drag numba → 0.66 and numpy → 2.4.x, which
  violate cudf/cuml/cugraph (`numba<0.61`) and cupy (`numpy<2.3`). The final pip step pins
  `numba>=0.59.1,<0.61` + `numpy>=2.0,<2.1` (→ numba 0.60.0 / numpy 2.0.2 / llvmlite 0.43.0).
  This drift did **not** itself cause the 0-GPU crash (the binding env var did), but running
  cuml/cugraph out of their supported numba range is a latent risk for the phenograph step,
  so keep it aligned. Must stay in the LAST pip step so it wins over scanpy's transitive pull.
- **segger's writer calls `pl.from_torch`, which polars<1.26 lacks (patched in via a
  shim).** After the GNN trains + predicts, `data/writer.py` builds the output frame with
  `pl.from_torch(tensor, schema=...)` (4 calls). polars only added `from_torch` **after 1.27**
  (PR pola-rs/polars#22177), but `cudf-polars 25.4` hard-caps polars at **`<1.26`** — so we
  can't get the real function without breaking the RAPIDS stack. Fix: keep polars at the
  RAPIDS-compatible **1.24** (pinned `>=1.24,<1.26`) and inject the missing function via a
  Docker `sed` step that appends, right after writer.py's `import torch`:
  `pl.from_torch = getattr(pl, "from_torch", lambda data, schema=None, **k:
  pl.from_numpy(data.detach().cpu().numpy(), schema=schema))`. The `from_numpy` equivalent was
  validated to reproduce segger's exact output (columns + Int64/Int64/Float32 dtypes) on the
  4-call + horizontal-concat pattern; `getattr` no-ops if a future polars ships a real one.
  The **path is hardcoded** (`/opt/conda/lib/python3.12/site-packages/...`, fixed by the
  py3.12 base) — a build failure there is a loud signal segger's layout moved. This surfaces
  only at the very end of a run (writer), so it was invisible until the GPU fix let training
  + prediction complete.
- **Empty (all-background) segmentation crashes `sd.to_polygons`.** If the incoming
  `segmentation.zarr` label image is entirely zero (no cells), `sd.to_polygons(...)` returns
  an **empty GeoDataFrame**, which then dies deep inside spatialdata's `ShapesModel.parse`
  (`data["geometry"].iloc[0]` → `IndexError: single positional indexer is out-of-bounds`) —
  before segger even runs. Root cause is upstream data, not segger: `2022_vizgen_human_lung_
  cancer_merfish_combined` has all-zero `cell_labels` from the stale **pre-0.8.0 spatialdata
  rasterize bug**, and `custom_segmentation` just copies that empty label (the same failure
  the pciseq method hit — see the `pciseq-empty-segmentation-lung-merscope` memory). Only this
  one dataset is affected (a prior survey of all processed datasets found no other empties), so
  it manifests as *only the `<this dataset> × custom_segmentation × segger` combo* failing.
  **Real fix: re-process that dataset under spatialdata≥0.8.0.** As defense, `script.py` now
  guards right after reading the segmentation (`int(seg_arr.data.max()) == 0 → raise` with an
  actionable message), mirroring pciseq's guard. Note the benchmark retries an exit-1 task 3×
  before ignoring it, so a guarded fail still burns 3 retries — re-processing is the cure.
- **Empty-`bd`-batch crash — a cell-free tile, NOT out-of-bounds transcripts (STILL OPEN).**
  segger builds a heterogeneous graph with `'tx'` (transcript) and `'bd'` (boundary/cell)
  node types and applies a `Positional2dEmbedder` to **both**. That embedder does
  `torch.zeros((batch.max()+1, 2))` (`ist_encoder.py`) — and `batch.max()` throws
  `RuntimeError: max(): Expected reduction dim to be specified for input.numel() == 0` on an
  **empty** tensor. A mini-batch drawn from **tiles that have transcripts but no boundaries**
  has zero `'bd'` nodes → empty `batch` → crash.
  - **Original hypothesis (now disproven as the cause):** the transcript-only tiles came from
    the combine step cropping labels to ~20000² while transcripts spanned a larger field
    (~38%, 7.6M/19.7M on mouse-brain rep3, sat outside the label). The adapter worked around it
    by EXCLUDING out-of-bounds transcripts from `transcripts.parquet`.
  - **What actually happened:** the crop bug was fixed at the source (`crop_points_by_global_xy`
    crops transcripts to the same global box as the labels). A full validation run (Seqera
    `1HlVGwXktVakGm`) then showed segger **completes on Xenium** (mouse brain rep3, human breast)
    but **still crashes on MERFISH** (`2022_vizgen_human_breast_cancer_merfish_combined/rep1`) —
    with `0/19222237 transcripts fall outside the label image` (**n_oob == 0**), training all 80
    batches, then dying at **validation tile ~21/28**. So the crash is **intra-field**: a tile
    covering tissue that has transcripts but no *segmented cells* (MERSCOPE fields have acellular
    / off-tissue gaps + a big extracellular transcript fraction; `custom_segmentation` copies the
    vendor mask, so transcript-only tiles exist). Xenium coverage is dense → never trips.
  - Because n_oob is 0, the OOB exclusion was a confirmed no-op and has been **removed**
    (replaced by the plain edge clamp; see step 3). It never addressed this crash anyway.
  - **The real fix (still to do — needs a container rebuild):** patch segger's
    `Positional2dEmbedder.forward` in `ist_encoder.py` to no-op on empty input, e.g. return
    `pos.new_zeros((pos.shape[0], 2 * self.dim))` when `pos.numel() == 0`, applied as a
    from_torch-style Docker `sed`/python shim in `config.vsh.yaml`. This is a **segger source
    patch (image rebuild)**, unlike the OOB simplification which is script-only. Verify the shim
    on a live `build_main` pod (`debug-component-k8s`) before baking it in — the GPU/RAPIDS
    rebuild is the expensive feedback loop. Not upstream in segger as of pinned `0233cf62`
    (PR #74 is a *different* fix — Shapely polygon validity, not this).
- **No CI test coverage** — GPU-only; can only be exercised on a CUDA host.
- **Small panels break segger's gene-embedding PCA.** segger's `setup_anndata` does
  `PCA(n_components=in_channels=128)` on a genes×genes correlation matrix built from genes
  that survive its count filters. On the test crop the 248-gene panel collapses to ~62
  surviving genes → `n_components=128 > 62` → `ValueError`. Handled by the `--node_
  representation_dim` arg + the PCA-ceiling auto-retry (step 4). This is a *small-data* limit,
  not a bug — full panels keep ≥128 genes. Watch: if a future run fails the PCA *after* the
  retry, segger's ceiling regex may have changed, or the surviving count is <2.
- **Depends on segger's CLI + output contract** (why segger is now commit-pinned):
  the `segger segment` flags (`--prediction-graph-buffer-ratio`, `--prediction-mode`,
  `--n-epochs`) and the output `segger_segmentation.parquet` columns (`row_index`,
  `segger_cell_id`, `segger_similarity`, `similarity_threshold` — **no `keep`** as of
  0.1.0). Any segger bump can change these; re-verify against the pinned commit.
- **Input is Xenium-layout** (`transcripts.parquet` + `{cell,nucleus}_boundaries.parquet`
  + an **`experiment.xenium`** marker). segger 0.1.0 **auto-detects the platform** from
  marker files (`_infer_platform` in `io/preprocessor.py`) — there is no `--platform` CLI
  flag. Its Xenium detector reads only `analysis_sw_version` from `experiment.xenium`, which
  selects the **v2+** loader (unassigned sentinel `"UNASSIGNED"`, what we write) vs the v1
  loader (`"-1"`). Field names otherwise match `XeniumTranscriptFields`/`XeniumBoundaryFields`.
- **The `experiment.xenium` version is propagated, not hardcoded.** The Xenium dataset
  loaders (`tenx_xenium`, `tenx_atera`) read the real `analysis_sw_version` from the raw
  `experiment.xenium` into `tables["table"].uns["xenium_analysis_sw_version"]`, which
  `process_dataset` carries through as `tables["metadata"].uns` of `raw_ist.zarr`.
  `script.py` reuses it when it is **≥2**; it falls back to a v2+ default (`xenium-3.0.0`)
  when the source is v1, missing, or **non-Xenium** (segger wraps MERFISH/CosMx as
  Xenium-layout too, so there's often no source version) — because we always write the
  `"UNASSIGNED"` sentinel, which only the v2+ loader understands. So the version is a
  fixed property of our *output* convention, informed by real provenance when available.
- **Pixel lookup truncates (not rounds)** transcript coords to int64 before the
  label-image index — matches `basic_transcript_assignment`'s convention.
- The Docker step ordering (esp. the final pandas pin) is fragile; changing the base
  image usually cascades into all the version pins above.

## Optimization / tuning

Parameter sweep for optimizing segger's transcript-assignment quality. The sweep files live
under `scripts/run_benchmark/param_sweep/` (`segger_params.yaml` + `run_test_segger_nebius.sh`).
`run_benchmark` expands the params file as a **star around the default** (one extra variant per
swept value, that ONE arg overridden — NOT a grid), so total variants =
1 default + Σ(sweep list lengths) = **8**. GPU runs are expensive, so the sweep is kept lean and
touches **only already-exposed args ⇒ SUBMITTABLE AS-IS (no container rebuild)**.

### Non-default audit (vs segger 0.1.0 @ `0233cf62`)

Every exposed quality knob's shipped default MATCHES segger 0.1.0's own default — there is **no**
pre-existing deviation to force onto a sweep axis:

| arg | shipped | segger 0.1.0 default | source |
|-----|---------|----------------------|--------|
| `n_epochs` | 20 | 20 | `cli/segment.py` (`n_epochs=20` literal) |
| `prediction_graph_buffer_ratio` | 0.05 | 0.05 | `data/data_module.py:151` (`= 0.05`) |
| `prediction_mode` | `cell` | `cell` | `data/data_module.py:149` (`prediction_graph_mode="cell"`) |
| `node_representation_dim` | 128 | 128 | `data/data_module.py:138` (`cells_embedding_size=128`) |

⚠️ **Correction to the arg-table note above:** the `0.5` mentioned there was the default of the
**OLD** `--prediction-expansion-ratio` flag. In 0.1.0 that flag was renamed to
`--prediction-graph-buffer-ratio` AND its default was lowered to **0.05** — exactly what we ship.
So 0.05 is *not* a deviation; it is segger's current default. (`transcripts_key` /
`coordinate_system` are our adapter's I/O keys, not segger knobs — fixed, never swept.)

### Tiers

- **Tier 0 — inputs, not parameters.** The biggest levers aren't in the sweep: the segmentation
  prior fed in (which segmentation method, how good its boundaries) and how densely those
  boundaries cover the transcript field (sparse coverage → cell-free tiles → the empty-`bd`
  gotcha). Both are fixed by the upstream stage here.
- **Tier 1 — highest impact on *what* gets assigned:**
  - **`prediction_mode`** {`nucleus`,`cell`,`uniform`} — which polygon set drives the prediction
    graph. We write the SAME polygons to both cell/nucleus files, so `nucleus` vs `cell` differ
    only in how segger treats them internally; `uniform` drops the prior-boundary identity.
    Default `cell`; sweep the other two → `[nucleus, uniform]`.
  - **`prediction_graph_buffer_ratio`** — fraction of each polygon's equivalent radius used to
    buffer (expand) it when building the prediction graph. The direct **recall↔precision** dial:
    a larger buffer expands each cell so it captures more surrounding transcripts (↑recall,
    ↓precision). Default 0.05 (tight); sweep UP → `[0.1, 0.25, 0.5]` (0.5 == the old
    expansion-ratio default — a natural upper anchor).
- **Tier 2 — quality/speed trade-off:**
  - **`n_epochs`** — GNN training epochs. More epochs = better convergence at a linear GPU-time
    cost (60 ≈ 3× the default's train time). Default 20; sweep UP → `[40, 60]`.
- **Tier 3 — not exposed (would need a config arg + a GPU/RAPIDS container rebuild — deliberately
  NOT done, to keep the sweep submittable-as-is).** Highest-value candidates for a serious
  follow-up: `prediction_max_k` (default 3 — k for the prediction kNN graph, another recall dial),
  `cells_representation` (`pca`/`morphology`), `segmentation_loss` (`triplet`/`bce`),
  `transcripts_max_dist`/`transcripts_max_k` (transcript-transcript graph radius), `learning_rate`.

### Why `node_representation_dim` is NOT swept

It maps to segger's `cells_embedding_size`/`in_channels` and sets the PCA `n_components` for the
gene embeddings. It is effectively an **upper bound**: on small panels/crops fewer genes survive
segger's count filters than 128, and `script.py` already auto-retries at the ceiling segger
reports. Sweeping it lower would only cap the embedding capacity (not improve quality); higher
would just re-hit the same ceiling. It is a *capacity/robustness* knob, not a quality axis → left
fixed at 128.

### Sweep summary

| tier | arg | default | swept values |
|------|-----|---------|--------------|
| 1 | `prediction_mode` | `cell` | `nucleus`, `uniform` |
| 1 | `prediction_graph_buffer_ratio` | 0.05 | 0.1, 0.25, 0.5 |
| 2 | `n_epochs` | 20 | 40, 60 |

Total = 1 default + 2 + 3 + 2 = **8 variants**. All args already exposed ⇒ **SUBMITTABLE AS-IS**.
segger is GPU-only: the `gpuh100` label in the config + `src/base/labels_nebius.config`
(`runAsUser:0` + Memory-backed `/dev/shm`) pin it to the GPU node group — no GPU-specific change
is needed in the sweep script.
