---
name: render-qmd-pod
description: >-
  Render a Quarto (.qmd) report against a given dataset-paths YAML on a pod on the
  kubernetes cluster, using a Viash component's image (default process_dataset) so the
  report's spatialdata/anndata/scanpy stack reads the combined datasets on the shared
  scratch PVC natively. Use when you need to render/execute a .qmd (e.g. the dataset
  visualisation report in scripts/results_report/) that reads data living ONLY on the
  cluster's tower-scratch PVC — a local render can't see those files. Handles pod
  creation, installing the render deps the image lacks (quarto, jupyter, spatialdata_plot),
  staging the qmd+yaml, rendering, and copying the self-contained HTML back.
---

# render-qmd-pod

The combined datasets (`/scratch/task_ist_preprocessing/datasets/…` +
`…/raw/…` single-cell refs) live **only** on the cluster's `tower-scratch` PVC, so a
report that reads them can't be rendered locally. This skill renders the `.qmd` on a
pod that **mounts that PVC** and runs a **component image** whose stack already reads
`.zarr`/`.h5ad` natively (default: `data_processors/process_dataset`, which ships
spatialdata 0.8 / anndata / scanpy / umap-learn / matplotlib). It then adds only the
report-rendering pieces the image lacks — **quarto**, **jupyter/ipykernel**, and any
report-specific python deps (**spatialdata_plot** by default) — stages the report,
renders, and copies the self-contained HTML back.

The report/YAML contract (as in `scripts/results_report/Dataset_report.qmd`): the qmd
reads its paths from `$DATASET_CONFIG` (falling back to `dataset_paths.yaml`), resolved
relative to the qmd's own directory, and the YAML's `resources_roots` points at the PVC
mount. `render` sets `$DATASET_CONFIG` to the staged YAML's basename and runs quarto with
cwd = the qmd's dir, so paths resolve exactly as they do locally.

## Prerequisites

- `kubectl` pointed at the cluster (`kubectl config current-context`). On Nebius the
  token expires — if a call prints "Switch to your browser to complete authentication",
  re-auth before the next long transfer/exec.
- A **dataset-paths YAML** whose `resources_roots` points at the PVC mount
  (`/scratch/task_ist_preprocessing`) and whose file paths exist there. To find the
  smallest self-contained dataset (full SC ref + panel-subset SC ref + spatial zarr all
  on the PVC), list `/scratch/task_ist_preprocessing/datasets` and `…/raw` on any pod
  that mounts `/scratch` (e.g. the long-lived `debug` pod). `dataset_paths_skin_k8s.yaml`
  is a worked example.
- The chosen component image must exist on ghcr (`check-component <comp>` if unsure).

## How to run

```bash
S=.claude/skills/render-qmd-pod/render_qmd.sh
$S up                                                   # create pod (scratch mounted) + install render deps
$S render <file.qmd> <paths.yaml> [out.html]            # stage both, render, copy HTML back (default ./<stem>.html)
$S down                                                 # delete the pod when done
$S sh '<bash>'                                          # poke the pod (debugging)
$S image                                                # print the resolved image ref
```

Example (the dataset report on the smallest self-contained dataset):

```bash
$S up
$S render scripts/results_report/Dataset_report.qmd \
         scripts/results_report/dataset_paths_skin_k8s.yaml \
         scripts/results_report/Dataset_report_skin_k8s.html
$S down
```

`up` is a one-time cost (~quarto 147MB download, cached on the PVC at
`/scratch/.render-cache`, + a pip install); after it, `render` reuses the pod and is
just the render itself. Re-render another dataset by pointing `render` at a different
YAML — no `up` needed. `up` is idempotent, so re-running it on an existing pod just
re-checks the deps.

## Adapting

- **Different image** (a report needing a different stack): `RENDER_COMPONENT=<ns/comp>`
  (resolved from `_viash.yaml` like the cluster pulls it) or `RENDER_IMAGE=<full ref>`.
- **Extra python deps** the image lacks: `RENDER_PIP="pkg_a pkg_b"` (default
  `spatialdata_plot`; the process_dataset image already has scanpy/spatialdata/umap-learn).
- **Other clusters**: `SCRATCH_PVC` / `SCRATCH_NS` / `SCRATCH_MOUNT` / `RENDER_POD` /
  `QUARTO_VERSION` / `RENDER_KERNEL`.

## Gotchas

- **The YAML must point at the PVC, not a local checkout.** Its `resources_roots` first
  entry must be the PVC mount (`/scratch/task_ist_preprocessing`); a local
  `~/projects/.../resources_test` path won't exist in the pod. Keep the local path only as
  a lower-priority fallback.
- **On-disk multiscale groups are `s0…s4`, not `scale0`.** spatialdata normalizes them to
  `scale0` keys in memory, so a report accessing `img["scale0"]["image"]` works unchanged —
  don't "fix" the report to `s0`.
- **`sleep infinity`, not a finite sleep** — the helper already does this; a completed pod
  can't be exec'd.
- **Scratch PVC is namespace-scoped** (`tower-nf`), so the pod runs there and all
  exec/cp/down use that namespace automatically.
- **If a render exec drops mid-stream**, the pod persists and the working copy + `.html`
  stay in `/tmp/render-qmd` inside the pod — re-run `render` (cheap once deps are installed)
  or fetch the file with `$S sh 'ls -la /tmp/render-qmd'` + `kubectl cp`.
- **Leaving the pod up** holds a PVC mount + a node slot. Fast to re-render, but run
  `$S down` when finished.

## After you render

- Open the HTML locally (`open <out.html>`) to check the figures.
- The self-contained HTML embeds all figures (base64 PNGs), so it's shareable as-is.
- `$S down` to release the pod.
