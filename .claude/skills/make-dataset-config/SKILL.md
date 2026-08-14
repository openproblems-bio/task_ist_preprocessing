---
name: make-dataset-config
description: >-
  Create a Dataset_report dataset-paths YAML (like results/config/*.yaml) for a combined
  dataset that lives on the cluster scratch PVC, and download its raw single-cell reference
  into raw/ if that .h5ad is missing on the pod. Use when you need a per-dataset config for
  the dataset visualisation report (the input to the render-qmd-pod skill), when asked to
  "make the config/yaml for <dataset> on the pod", or to stage a dataset's raw SC reference.
  Writes results/config/<id>.yaml with an id matching the dataset's id on the pod exactly,
  auto-resolving the raw SC source from the combine scripts and pulling it from the public
  openproblems-data S3 bucket when absent.
---

# make-dataset-config

The dataset visualisation report needs three inputs — the **full** single-cell reference,
the **panel-subset** SC reference, and the **spatial** zarr. On the cluster the panel-subset
SC (`output_sc.h5ad`) and spatial (`output_sp.zarr`) already sit in the combined dataset dir,
but the **full raw SC reference** (the combine step's `input_sc`) is often **not staged** on
the PVC — only custom references (ltx/mpii/kuppe/ganier) are. This skill fills that gap: it
resolves the raw SC path from the combine scripts, **downloads it into `raw/` if missing**,
and writes a ready-to-render `results/config/<id>.yaml`. Pair it with **render-qmd-pod** to
produce the HTML.

Why the download is safe & 1:1: standard datasets set `input_sc` to
`$input_dir/<subpath>` with `$input_dir=s3://openproblems-data/resources/datasets` — a public
bucket — so `<subpath>` maps exactly onto the PVC (`raw/<subpath>`) and onto the HTTPS endpoint
(`https://openproblems-data.s3.amazonaws.com/resources/datasets/<subpath>`).

## Prerequisites

- `kubectl` on the cluster (`kubectl config current-context`). On Nebius the token expires —
  re-auth if a call prints "Switch to your browser to complete authentication".
- A pod that **mounts the scratch PVC** at `/scratch` and has `curl` + `python3`/`h5py`
  (default `debug` in `tower-nf`). Override with `SCRATCH_POD`; if none exists, bring one up
  via the render-qmd-pod skill (`$S up`) and set `SCRATCH_POD=render-qmd`.

## How to run

```bash
S=.claude/skills/make-dataset-config/make_dataset_config.sh
$S <dataset_id>                              # auto-resolve raw SC, download if missing, write config
$S <dataset_id> --sc-source <subpath>        # override the S3 subpath under resources/datasets/
$S <dataset_id> --sc-raw <raw-rel-path>      # raw SC already on the PVC (custom ref); no download
$S <dataset_id> --out <file>                 # custom output path
```

`<dataset_id>` is the dataset's dir under `datasets/` on the pod — the id in
`results/raw_dataset/datasets.csv` (e.g. `2024_10x_human_liver_xenium_combined`, or a rep like
`2023_10x_mouse_brain_xenium_combined/rep1`). Find ids with
`kubectl exec -n tower-nf debug -- ls /scratch/task_ist_preprocessing/datasets`.

Example (the healthy-liver dataset — auto-resolves `2022Andrews_human_liver_sc`, downloads it
to `raw/` if absent, writes `results/config/2024_10x_human_liver_xenium_combined.yaml`):

```bash
$S 2024_10x_human_liver_xenium_combined
```

Then render it:

```bash
R=.claude/skills/render-qmd-pod/render_qmd.sh
$R up
$R render results/Dataset_report.qmd results/config/2024_10x_human_liver_xenium_combined.yaml \
          results/raw_dataset/2024_10x_human_liver_xenium_combined.html
$R down
```

## What it does

1. Verifies `datasets/<id>/{output_sp.zarr,output_sc.h5ad}` exist on the pod.
2. Resolves the raw SC reference: parses `input_sc` for that id from
   `scripts/create_resources/combine/*.sh`; `$input_dir/<subpath>` → `raw/<subpath>` + S3 URL.
   (If it can't — e.g. a `$sc_ref` variable — it stops and asks for `--sc-source`/`--sc-raw`.)
3. Downloads the `.h5ad` into `raw/` **only if missing** (idempotent; atomic `.part`→final).
4. Best-effort checks the ref has `X_pca` + `normalized` + `cell_type` (what the report reads).
5. Writes `results/config/<id>.yaml` (slashes in reps → `_`), with `id:` = the exact pod id,
   `scrnaseq_raw` = the raw ref, and `dataset.dir` = the combined dir.

## Gotchas

- **The id must match the pod dir exactly.** Rep datasets carry a `/repN` (or `/B_…`,
  `/978_reg1`) suffix — pass the full leaf id. The output filename slugifies `/`→`_`, but the
  in-file `id:` and `dataset.dir` keep the true path.
- **Custom references (ltx/mpii/kuppe) aren't on S3** — their `input_sc` is a pod-local/`$sc_ref`
  path, so auto-resolve can't form a URL. They're usually already on the PVC: pass
  `--sc-raw raw/<their-path>/dataset.h5ad` (find it with
  `kubectl exec -n tower-nf debug -- find /scratch/task_ist_preprocessing/raw -name dataset.h5ad`).
- **Big downloads land on shared storage.** Refs range ~380 MB–1.5 GB; the download is skipped
  if the file is already there, so re-running is cheap.
- **`raw/` mirrors the S3 layout** (`raw/scrnaseq_for_ist/<name>/dataset.h5ad`) so a later
  pipeline run or another config reuses the same file — don't invent a different raw path.

## After you create the config

- Render it with **render-qmd-pod** to confirm the ref + paths resolve end-to-end.
- The config is a repo file under `results/config/`; the downloaded `.h5ad` stays on the PVC.
