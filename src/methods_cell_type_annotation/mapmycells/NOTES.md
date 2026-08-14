# mapmycells — implementation & tuning notes

Cell-type annotation by mapping a single-cell reference onto the spatial data with the
**Allen Institute `cell_type_mapper`** (the engine behind the online **MapMyCells** tool).
Despite the config `summary` ("...using moscot"), this component does **not** use moscot —
that line is a copy-paste leftover from the moscot component. It shells out to the Allen
`cell_type_mapper` CLI.

## Citation caveat

`references.doi` = `10.1038/s41586-023-06812-z` → Yao et al. 2023, *Nature*, "A high-resolution
transcriptomic and spatial atlas of cell types in the whole mouse brain" — the atlas paper that
introduced/used the MapMyCells hierarchical-mapping algorithm. This is a reasonable citation for
the *algorithm*; the *software* is documented in the repo (no standalone method paper). Fine to
keep; optionally also cite the repo.

## How it works (script.py)

Two `cell_type_mapper` CLI calls, shelled out via `subprocess.run`:

1. **`cli.precompute_stats_scrattch`** — builds reference statistics from the SC reference.
   - `--h5ad_path` = SC reference (with `adata_sc.X` set to the `counts` layer if present).
   - `--hierarchy "['cell_type']"` — a **single-level (flat) taxonomy** keyed on the reference's
     `cell_type` column.
2. **`cli.from_specified_markers`** — the actual correlation mapping of the query (spatial) cells.
   - `--query_path` = spatial `input_spatial_normalized_counts`.
   - `--type_assignment.normalization log2CPM` — query is log2CPM-normalized before correlating.
   - `--query_markers.serialized_lookup` = `genes.json` = `{"None": common_genes}` — i.e. the
     **entire intersected gene panel is used as markers** (no marker selection). This is Tier-0
     input behaviour, not a knob: spatial panels are small, so all shared genes are "markers".
   - `--flatten True` — flatten the (already single-level) taxonomy.
   - `--type_assignment.bootstrap_iteration` / `--type_assignment.bootstrap_factor` — bootstrap
     controls (see below). **These were hardcoded `1` and `1.0`; now exposed as config args.**

Output: `results.csv`'s `cell_type_label` column is written to `adata_sp.obs[celltype_key]`.

### Input contract / genes
Only genes present in **both** the spatial panel and the SC reference are used
(`common_genes`), for both the precompute and the query. Reference counts come from the
`counts` layer when present, else `.X`.

## Environment / build

- Image: `openproblems/base_python:1` + `setup_spatialdata_partial` + `setup_txsim_partial`,
  plus pip `numpy` and `git+https://github.com/AllenInstitute/cell_type_mapper.git` (installed
  from `main`, **unpinned** — upstream CLI/param changes could drift).
- **CPU-only.** Nextflow directives `label: [ hightime, midcpu, highmem ]` — no `gpu` label,
  no torch/CUDA in the script. It's a statistical correlation mapper.

## Optimization / tuning

The mapper assigns each query cell by correlating it against the reference cell-type profiles.
**Bootstrapping** repeats that assignment `bootstrap_iteration` times, each on a random
`bootstrap_factor` fraction of the marker genes, and majority-votes the label — this is the
robustness mechanism (and produces a per-cell confidence fraction).

Documented settings (Allen `docs/mapping_cells.md`):
- **No bootstrapping (fastest):** `bootstrap_iteration 1`, `bootstrap_factor 1.0` — "an order
  of magnitude faster".
- **Bootstrapped correlation mapping (robust):** `bootstrap_iteration 100`, `bootstrap_factor 0.9`.

**The component shipped the fast, no-bootstrap setting (1, 1.0)** — a deliberate speed choice.
Tuning here = walking `bootstrap_iteration` back up toward the tool default of 100.

### Tiers
- **Tier 0 (input, not a knob):** all common genes used as markers; single-level `cell_type`
  taxonomy; `log2CPM` query normalization (correctness — must match the reference space; left
  hardcoded, not swept).
- **Tier 1 — `bootstrap_iteration` (newly exposed).** The recall/robustness↔speed dial. Shipped
  `1` (no bootstrap) vs tool default `100`. **Swept: [10, 25, 50, 100].**
- **Tier 1/2 — `bootstrap_factor` (newly exposed).** Marker-downsampling fraction per iteration.
  Shipped `1.0` (no downsampling) vs tool default `0.9`. **Coupled:** with `factor = 1.0` every
  bootstrap iteration is identical, so the iteration axis would be inert. It is therefore **not
  given its own sweep axis**; instead the sweep's `default:` block pins it to the tool default
  `0.9` so the `bootstrap_iteration` sweep is meaningful. (A joint iteration×factor sweep is
  possible future work but violates the benchmark's one-arg-at-a-time star expansion.)
- **Performance only (fixed, never swept):** `n_processors`, `chunk_size`, `rng_seed`.

### Non-default audit
- `bootstrap_iteration`: config default `1` ≠ tool default `100` → on the sweep axis
  (`[10, 25, 50, 100]`, straddling shipped 1 → tool 100). ✔
- `bootstrap_factor`: config default `1.0` ≠ tool default `0.9` → the sweep's `default:` block
  sets it to `0.9` (the tool default), removing the deviation, because it can only take effect
  jointly with `bootstrap_iteration > 1` (see Tier 1/2). ✔
- Config defaults are deliberately left at the **shipped (1, 1.0)** so a full-benchmark run that
  does not pass these args reproduces the existing mapmycells baseline.

### Sweep
`scripts/run_benchmark/param_sweep/mapmycells_params.yaml`:
- `default: { bootstrap_iteration: 1, bootstrap_factor: 0.9 }`  (fast anchor, but with
  downsampling on so the axis is live)
- `sweep: { bootstrap_iteration: [10, 25, 50, 100] }`
- **Total = 1 default + 4 = 5 variants.**

⚠️ These args are **newly added**, so the `build/main` container must be **rebuilt + pushed**
before the sweep run understands them (`viash ns build` → container rebuild → push;
see the `check-component` skill). NOT submittable until the image is refreshed.

Runtime: `bootstrap_iteration 100` is the slowest variant (~10× the no-bootstrap pass per the
docs); on the small test datasets this is expected to be minutes, but it is the cost driver.

## Iteration log
- 2026-08-12 — Exposed `bootstrap_iteration` + `bootstrap_factor` (were hardcoded 1/1.0);
  wrote this NOTES.md; created the 5-variant sweep (iteration ladder, factor pinned 0.9).
  Config `viash config view` OK, `script.py` parses. Not yet built/pushed/run.
