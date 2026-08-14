# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an **OpenProblems v2 benchmarking task** for preprocessing imaging-based spatial transcriptomics (iST) data. It uses **Viash** (containerized component framework) and **Nextflow** (workflow orchestration) to benchmark methods across a multi-stage preprocessing pipeline: segmentation → transcript assignment → count aggregation → QC filtering → cell volume calculation → normalization → cell type annotation → expression correction.

## Key Commands

### Build & Test

```bash
# Build all components (creates executables in target/)
viash ns build --parallel

# Build Docker containers (required before running benchmarks)
scripts/project/build_all_docker_containers.sh
# Equivalent: viash ns build --parallel --setup cachedbuild

# Test a single component
viash test src/methods_segmentation/binning/config.vsh.yaml

# Test all components in parallel
viash ns test --parallel
# or via script:
scripts/project/test_all_components.sh

# Filter tests by namespace
viash ns test --parallel -q methods_segmentation
```

### Run Benchmarks

```bash
# Local test run (uses test resources, outputs to temp/results/)
scripts/run_benchmark/run_test_local.sh

# Full local benchmark (uses resources/, outputs to resources/results/)
scripts/run_benchmark/run_full_local.sh

# Cloud runs (Nebius, Seqera)
scripts/run_benchmark/run_test_nebius.sh
scripts/run_benchmark/run_test_seqeracloud.sh
```

### Create Test Resources

```bash
# Download and process test datasets (run via Nextflow)
scripts/create_test_resources/<dataset_script>.sh
```

## Architecture

### Component Structure

Every method/metric is a self-contained Viash component:

```
src/<namespace>/<component_name>/
  config.vsh.yaml    # Component spec: arguments, Docker image, resources
  script.py          # Implementation (or script.R)
```

`config.vsh.yaml` merges the API interface via `__merge__: /src/api/comp_<type>.yaml`, which enforces standard input/output arguments. Component-specific parameters go in the `arguments:` block.

### API Layer (`src/api/`)

- `comp_method_*.yaml` — Defines the standard interface (input/output args) for each pipeline stage
- `file_*.yaml` — Defines file format specifications
- Components inherit these via `__merge__` in their config

### Data Formats

All spatial data flows as `.zarr` (SpatialData format); single-cell references and outputs use `.h5ad` (AnnData). Key file types:
- `raw_ist.zarr` → raw imaging data
- `segmentation.zarr` → cell label maps
- `transcript_assignments.zarr` → per-transcript cell assignments
- `spatial_aggregated_counts.h5ad` → cell × gene matrix with centroids/volume
- `scrnaseq_reference.h5ad` → reference atlas for annotation

### Pipeline Flow

The benchmark workflow (`src/workflows/run_benchmark/main.nf`) fans out across all enabled methods at each stage, running every combination. The dataset processor (`src/workflows/process_datasets/main.nf`) converts raw vendor outputs to the standardized `raw_ist.zarr` format.

### Key Dependencies

- `txsim` (from `theislab/txsim@dev`) — core spatial transcriptomics utilities used in most scripts
- `spatialdata>=0.7.3` — standardized spatial data container
- `anndata>=0.12.0`, `zarr>=3.0.0`

Base Docker setups are mixed in via `__merge__` in component configs:
- `src/base/setup_spatialdata_partial.yaml` — spatialdata + dependencies
- `src/base/setup_txsim_partial.yaml` — txsim library

### Resource Labels (Nextflow)

Components are tagged with resource labels in their `runners[nextflow].directives`:
- Memory: `lowmem` (20 GB), `midmem` (50 GB), `highmem` (100 GB)
- CPU: `lowcpu` (5), `midcpu` (15), `highcpu` (30)
- Time: `lowtime` (1h), `midtime` (4h), `hightime` (8h), `veryhightime` (24h)

### Test Resources

Located in `resources_test/` — downloaded from S3 (`s3://openproblems-data/resources_test/`). The `_viash.yaml` defines the S3 paths; use `viash run` with the appropriate config to sync them.

### Common Submodule

`common/` is a git submodule from `openproblems-bio/openproblems`. It contains shared component test scripts (`common/component_tests/`) that validate output format compliance.

## Adding a New Method

1. Create `src/<namespace>/<method_name>/config.vsh.yaml` — merge the appropriate `comp_method_*.yaml` API
2. Create `script.py` implementing the logic; read inputs via `par["input_*"]`, write outputs via `par["output_*"]`
3. Add `test_resources` pointing to `resources_test/` data and a `common/component_tests/` test runner
4. Register the method in `scripts/run_benchmark/config.yaml` and the run scripts
