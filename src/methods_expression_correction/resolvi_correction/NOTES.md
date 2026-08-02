# resolvi_correction — NOTES

Authoritative how-it-works / why-the-setup / where-it-breaks reference for the
`resolvi_correction` expression-correction method. Committed with the code.

## What this component is

- **Stage / API:** `methods_expression_correction` (merges
  `/src/api/comp_method_expression_correction.yaml`). Runs late in the iST pipeline: it takes
  `spatial_with_cell_types.h5ad` (a cell x gene matrix that already has cell-type labels,
  centroids, and a `normalized` layer) and rewrites its `counts` / `normalized` layers with a
  denoised, spillover-/background-corrected version.
- **Underlying tool:** **resolVI** (`scvi.external.RESOLVI`) from scvi-tools — a variational
  autoencoder that models each observed cell's expression as a mix of (α0) true self-expression,
  (α1) diffusion/spillover from spatial neighbours, and (α2) unspecific background, then samples
  the posterior to reassign misplaced molecules.
- **GPU:** deep VAE model → tagged `gpu` in the nextflow directives. NB the docker `image:` is
  currently `openproblems/base_python:1` (CPU base) with an explicit TODO to move to a GPU
  pytorch image (`base_pytorch_nvidia:1`) — see the commented line in `config.vsh.yaml`. So the
  component requests a GPU node but the container is a CPU-first python image with `scvi-tools`
  pip-installed; whether it actually lights up CUDA depends on the torch build in that image.
- **Links:** docs https://docs.scvi-tools.org/en/latest/user_guide/models/resolvi.html ,
  repo https://github.com/scverse/scvi-tools .

## Citation

`references.doi: 10.1101/2025.01.20.634005` is **correct** (verified via bioRxiv): Ergen & Yosef,
2025, *"ResolVI - addressing noise and bias in spatial transcriptomics"* — the method-specific
paper, not a generic framework citation. No caveat needed.

## script.py — step by step

1. Read `input_spatial_with_cell_types`; stash the incoming `normalized` layer as
   `normalized_uncorrected` (script.py L32-33).
2. `sc.pp.filter_cells(min_genes=5)` — drop near-empty cells (L36).
3. Build the spatial neighbour input: stack `obs.centroid_x/centroid_y` into
   `obsm["X_spatial"]` (L38-39) — this is what resolVI uses to find neighbours.
4. `RESOLVI.setup_anndata(labels_key=celltype_key, layer="counts")` (L45).
5. Construct `RESOLVI(..., semisupervised=True, n_hidden=, encode_covariates=, downsample_counts=)`
   (L48-51). `semisupervised=True` is **hard-coded** (tutorials recommend the supervised variant
   when cell-type labels are available, which they are at this stage).
6. `train(max_epochs=100)` — **hard-coded** (L52).
7. Two `sample_posterior` calls (**`num_samples=20`, `batch_size=4000` both hard-coded**):
   - `model_corrected` → `px_rate` summarised at the q50/median (L55-59).
   - `model_residuals` → `mixture_proportions, mean_poisson, per_gene_background,
     per_neighbor_diffusion, px_r_inv` (L63-70). The exact `return_sites` set matters — see the
     scvi-tools issue #3252 link in the code.
8. Corrected counts = `counts * px_rate_q50 / (1 + px_rate_q50 + mean_poisson)` (L80-81), then a
   size-factor renormalisation into the `normalized` layer (L84-86). The renormalisation is
   acknowledged as non-ideal in a long code comment (L88-97): the pipeline runs normalization
   *before* correction, so the corrected counts are re-normalised here by a crude size factor
   rather than re-running the real normalization method.

**Input-contract note:** the code comment at L28 claims `par['input_sc']` is required, but the
script never reads `input_scrnaseq_reference` — only `input_spatial_with_cell_types`. The SC
reference is optional in the API and unused here.

## Arguments (exposed in config.vsh.yaml)

| arg | type | config default | upstream (scvi RESOLVI) default | notes |
|-----|------|----------------|--------------------------------|-------|
| `--celltype_key` | string | `cell_type` | n/a | `labels_key` for `setup_anndata`; input mapping, not a tuning knob (Tier 0) |
| `--n_hidden` | integer | `32` | **32** (matches) | VAE hidden-layer width; model capacity |
| `--encode_covariates` | boolean | `false` | `False` (matches; via `**model_kwargs`) | feed batch/covariates through the encoder |
| `--downsample_counts` | boolean | `true` | `True` (matches) | subsample per-cell counts to equalise depth in training |

**Hard-coded (NOT exposed) knobs that matter:** `max_epochs=100`, `num_samples=20` (upstream
`sample_posterior` default is **1000**), `mixture_k` (upstream default **50**, unset here),
`batch_size=4000`, `semisupervised=True`, `n_latent=10`, `n_layers=2`, `dropout_rate=0.05`.

## Optimization / tuning

Verified upstream defaults (scvi-tools RESOLVI): `n_hidden=32`, `n_hidden_encoder=128`,
`n_latent=10`, `n_layers=2`, `dropout_rate=0.05`, `mixture_k=50`, `downsample_counts=True`;
`sample_posterior(num_samples=1000)`.

**Non-default audit result:** all three exposed *tunable* args (`n_hidden`, `encode_covariates`,
`downsample_counts`) sit at their scvi upstream defaults — nothing was silently deviated, so
nothing is *forced* onto a sweep axis. They are swept because they are the only exposed levers.

Tiers:

- **Tier 0 — input, not a parameter.** `celltype_key` (which obs column carries the labels),
  the spatial neighbour graph built from `centroid_x/y`, and the `semisupervised=True` choice.
  Biggest lever overall is segmentation/assignment quality upstream, not resolVI's own knobs.
- **Tier 1 — highest posterior-quality impact, but NOT EXPOSED (see below).** `num_samples`
  (posterior samples; script uses **20** vs upstream **1000** — the single biggest
  quality-vs-time dial for denoising stability) and `max_epochs` (**100**, hard-coded).
- **Tier 2 — exposed quality/capacity trade-offs (what the sweep varies).**
  - `n_hidden`: 32 (default) → **64, 128**. More capacity to model diffusion + background,
    slower. (script.py's own grid-search note lists exactly 32/64/128.)
  - `encode_covariates`: false → **true** (boolean flip).
  - `downsample_counts`: true → **false** (boolean flip).
- **Tier 3 — high-value knobs NOT surfaced by the component (future work; would need a config
  arg + container rebuild — NOT in the submittable sweep):**
  - `num_samples` (posterior samples): expose as e.g. `--num_samples`, sweep `[50, 100, 200]`
    straddling the script's 20 toward upstream's 1000. Highest expected quality gain.
  - `max_epochs`: expose as `--n_epochs`, sweep around 100 (e.g. `[50, 100, 200, 400]`).
  - `mixture_k`: number of background/diffusion mixture components (upstream 50); expose and
    sweep e.g. `[10, 50, 100]`.

  Exposing any of these needs: add to `config.vsh.yaml arguments:`, thread `par[...]` into the
  matching `RESOLVI(...)` / `train(...)` / `sample_posterior(...)` call in `script.py`, then
  `viash ns build` + a container rebuild (see `check-component`) before a build/main Nebius run
  would honour them.

## Wiring

- Enabled at the `expression_correction` stage. Present (commented) in the stock run scripts
  (`run_test_nebius.sh`), active in `run_gpu_nebius.sh` alongside `no_correction`.
- Sweep: `scripts/run_benchmark/param_sweep/resolvi_correction_params.yaml` (committed params)
  + `scripts/run_benchmark/param_sweep/run_test_resolvi_correction_nebius.sh` (Nebius launcher,
  GPU compute env `5hfmdCBxMRd4nHZaJKYEQZ`, labels include `gpu`). Total variants =
  1 default + 2 (`n_hidden`) + 1 (`encode_covariates`) + 1 (`downsample_counts`) = **5**.

## Risk points / gotchas

- **Submittable-now constraint:** the Nebius run uses `--revision build/main`, so the sweep can
  only vary args already in the build/main container (the four above). Any newly-exposed Tier-3
  arg is silently ignored until rebuilt.
- **`return_sites` set is load-bearing** (script.py L63-70, scvi issue #3252) — changing it can
  break the residual sampling.
- **CPU image / gpu label mismatch** — see "What this component is". Confirm the image actually
  provides a CUDA torch before relying on GPU acceleration.
- **Renormalisation is crude** (script.py L84-97) — the corrected `normalized` layer is a
  size-factor rescale, not a re-run of the pipeline's real normalization step.
