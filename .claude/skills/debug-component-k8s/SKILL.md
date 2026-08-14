---
name: debug-component-k8s
description: >-
  Troubleshoot a Viash component/module's RUNTIME environment by running its container
  image live on the kubernetes cluster — dependency conflicts, missing/mismatched packages,
  native-library (.so / UCX / cv2) errors, CLI/API changes in a pinned tool, base-image
  problems. Use when a component fails INSIDE its image at run time (ImportError,
  ModuleNotFoundError, version TypeError, undefined symbol, SIGABRT, a CLI "unknown option"),
  and you need to reproduce and fix it without a full Nextflow/Seqera run. Complements
  check-component (which only checks whether the image is deploy-fresh, not whether it works).
---

# debug-component-k8s

`check-component` tells you *whether the right image is deployed*. This skill tells you
*whether that image actually works* — by running it on the cluster and poking at its
environment directly, so you find and validate a fix **before** paying for a ~30-min
rebuild and a GPU benchmark run.

The core loop: **spin up a throwaway pod on the component's image → reproduce the failure
with small probes → iterate on a fix live (pip/conda/env) → encode the validated recipe in
`config.vsh.yaml` → rebuild.** Each real benchmark run and each rebuild is expensive, so
every fix should be proven in the pod first.

## Prerequisites

- `kubectl` pointed at the cluster (`kubectl config current-context`). On Nebius the token
  expires — if a call prints "Switch to your browser to complete authentication", re-auth
  (`nebius ...`) or the next long transfer/exec dies mid-stream.
- The component's image must exist on the registry (run `check-component <comp>` first; if
  it's stale/missing, fix that before debugging the wrong image).

## How to run

```bash
S=.claude/skills/debug-component-k8s/inspect_pod.sh
$S up   <component>            # launch pod on <comp>:build_main, wait Ready (no GPU)
$S up   <component> --gpu      # ...on a GPU node (only for rmm/cuInit or a real run)
$S up   <component> --scratch  # ...with the shared scratch PVC mounted at /scratch
$S run  <component> '<bash>'   # run a bash snippet inside the pod
$S py   <component> '<python>' # run python (RAPIDS_NO_INITIALIZE preset so cudf imports)
$S down <component>            # delete the pod when done
```

Most dependency / import / native-lib / CLI troubleshooting needs **no GPU** — the default
CPU pod schedules anywhere and is fast. Only reach for `--gpu` to exercise code that touches
the device.

## Inspecting a run's data files (`--scratch`)

When you need to look at the actual `.h5ad` / `.zarr` a task staged or produced (shapes,
`.X`, `obs`/`var`, layers, why an assertion tripped), **run the component's OWN image** —
its stack (anndata, spatialdata, txsim, …) reads those files natively, so you never
pip-install readers into a generic pod. `--scratch` mounts the shared scratch PVC
(`tower-scratch`) at `/scratch`; because PVCs are namespace-scoped the pod then runs in the
PVC's namespace (`tower-nf`), and `run`/`py`/`down` find it there automatically.

```bash
$S up moscot --scratch
$S py moscot 'import anndata as ad; a = ad.read_h5ad("/scratch/<hash>/_viash_par/input_scrnaseq_reference_1/output_sc.h5ad"); print(a.X.shape, a.X.dtype, list(a.layers))'
$S down moscot
```

- **You must know the exact work-dir path — ASK THE USER for it.** Paths are opaque Nextflow
  hashes (`/scratch/<runhash>/<xx>/<taskhash>`); you cannot guess them. Given one, list it
  first to orient (via the running pod, or a live `debug` pod that already mounts `/scratch`).
- Inside a task work dir, **inputs are staged under `_viash_par/<arg_name>_1/…`** (one dir per
  API input; the file keeps its upstream name) and the **output** under the published
  `mpii_.../cta_<method>.<method>.output.h5ad` subpath. `.command.sh` lists the exact
  `VIASH_PAR_*` paths, and `.command.err` / `.exitcode` show how the task ended.
- Pick the image whose deps match the file: an `.h5ad` → any component image (all have
  anndata); a SpatialData `.zarr` → a component that installs `spatialdata`.
- For other clusters, override `SCRATCH_PVC` / `SCRATCH_NS` / `SCRATCH_MOUNT` in the env.
- Reading is fine even on non-root images — staged files are world-readable.

## The diagnostic playbook

Reproduce the actual error first, then work outward. Useful probes (adapt the package/module
names to the failure):

- **Reproduce the import/run:** `$S py <comp> 'from X import Y'` — get the FULL traceback
  (don't `2>/dev/null` it away). The deepest frame names the culprit module/symbol.
- **Versions + who-requires-what:** compare a package's installed version to what its
  siblings *demand* — this is how you spot version conflicts:
  ```python
  import importlib.metadata as m
  print(m.version("pandas"))
  print(m.requires("anndata"))   # -> ['pandas>=2.3', ...]  (the conflicting pin)
  ```
- **Is it installed at all?** `python -c "import cugraph"` per module → `MISSING`/`OK`.
- **Read the pinned tool's source** (contracts drift on unpinned installs): grep the
  installed package for the CLI/schema you depend on, e.g.
  `grep -rn "segger_cell_id\|write_parquet" $(python -c "import segger,os;print(os.path.dirname(segger.__file__))")`.
  Check a CLI's current options with `<tool> <subcmd> --help`.
- **Native-library errors** (`libGL.so.1: cannot open`, `undefined symbol`, `double free`):
  `ldd <lib>.so`, `nm -D <lib>.so | grep <symbol>`, `find / -name 'libX.so*'` to spot a
  **system lib shadowing a wheel's lib** (the classic cause — e.g. `/opt/hpcx/ucx` vs a pip
  UCX wheel). Fix via the missing system pkg, a headless variant of the wheel, or
  `LD_LIBRARY_PATH` ordering — but if it's a *heap* corruption (SIGABRT/double-free), a pin
  won't save you; the base image is wrong (see below).
- **Validate a candidate fix live** before editing the config: `$S run <comp> 'pip install
  "anndata>=0.12,<0.13" ... && python -c "..."'`. Re-check that the fix didn't disturb other
  pins (`importlib.metadata.version` before/after).

## Gotchas (learned the hard way)

- **`sleep infinity`, not `sleep N`.** A finite sleep (or the image's own entrypoint) can let
  the pod reach phase `Succeeded`; you then get `cannot exec into a completed pod`. The
  helper's `up` already uses `sleep infinity`.
- **Long-lived `kubectl exec` streams time out** (`i/o timeout` / "reading from error stream")
  over a flaky link. Keep execs SHORT. For long work (a build, a `docker save`, a poll),
  detach it inside the pod with `setsid`/`nohup ... < /dev/null &` and write an exit-code
  marker (`; echo $? > /work/done`), then poll the marker with short execs. `pgrep -f "cmd"`
  false-matches your own poll command — check for a done-marker/file-size instead.
- **CPU-only pod ⇒ GPU init fails, and that's EXPECTED.** `numba: CUDA driver library cannot
  be found`, `Failed to dlopen libcuda.so.1`, or a traceback ending in `rmm.reinitialize` /
  `cuInit` just means *no GPU on this pod* — not a real bug. Distinguish it from genuine
  errors; set `RAPIDS_NO_INITIALIZE=1 CUDF_NO_INITIALIZE=1 RMM_NO_INITIALIZE=1` (the `py`
  helper does) so `import cudf` and friends don't abort. To test past it, use `--gpu`.
- **Non-root base images bite twice — at build AND at run.** Check `whoami`/`id -u`.
  `rapidsai/base` runs as `rapids` (uid 1001).
  - *Build:* a viash `type: apt` step FAILS (needs root). Install system-ish deps via conda
    (`conda install -c conda-forge git`) and prefer wheels that avoid system libs (e.g.
    `opencv-python-headless` instead of `opencv-python`, which drops the libGL dependency).
  - *Run:* the container can't write Nextflow's **root-owned task scratch dir** →
    `Permission denied` on `.command.log` / `.command.begin` / `.exitcode` (fails before your
    script runs). viash can't set the image `USER`, so force root at the Nextflow level: add
    `[runAsUser: 0]` to that process's `pod` directive (e.g. on the `gpu` label in
    `labels_nebius.config`). Harmless for already-root images.
- **Don't infer from timestamps/labels — get ground truth.** When two signals disagree (a
  label says current, the run says stale), *run the actual image* and read versions/imports.
  Registry `created` timestamps are UTC; git commit times carry a zone — comparing them
  raw is a classic false conclusion.
- **Pin unpinned upstream installs.** A `github: owner/repo` (no `@commit`) install tracks
  trunk; the tool's CLI/output schema can change under your `script.py`. Pin to the commit you
  validated (`owner/repo@<sha>`; get it from
  `cat $(...)/segger-*.dist-info/direct_url.json`).
- **Heap corruption / SIGABRT from mixed native stacks ⇒ change the base image, not the pins.**
  If a required library stack (e.g. full RAPIDS) fundamentally conflicts with the base's own
  native libs (e.g. NGC's HPC-X UCX), stop patching — rebase on an image built for that stack
  (e.g. `rapidsai/base`) and layer the rest on top. Validate the whole layered stack imports
  (exit 0) in the pod before committing the new base.

## After you find the fix

1. Encode the validated recipe (pins, base image, extra installs) in the component's
   `config.vsh.yaml` setup steps — **last-wins ordering matters** for pins that must beat a
   github install's unpinned deps.
2. `viash config view src/<ns>/<comp>/config.vsh.yaml >/dev/null` to confirm it parses.
3. `$S down <component>` to clean up the pod.
4. Commit → merge to `main` → CI rebuilds → `check-component <comp>` confirms the new image is
   live → re-run. The pod only proves imports/coexistence; the GPU benchmark run is the final
   confirmation for anything device-dependent.
