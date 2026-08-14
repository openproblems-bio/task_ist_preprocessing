# iST Preprocessing Benchmark — Pipeline Elements

> Reference sheet describing every element of the `task_ist_preprocessing` benchmark, for use as an input to visualization/design work. Organized so each block (nodes, stages, methods, formats, skills) can be turned directly into diagram elements.

---

## 0. What this is

An **OpenProblems v2 benchmarking task** for preprocessing **imaging-based spatial transcriptomics (iST)** data (10x Xenium, Vizgen MERSCOPE, Bruker CosMx, Allen MERFISH).

**Tech stack (the "rails" everything runs on):**

| Layer | Tool | Role |
|-------|------|------|
| Component framework | **Viash** | Each method/metric is a containerized component (`config.vsh.yaml` + `script.py`/`.R`) |
| Orchestration | **Nextflow** | Fans components out across the benchmark matrix |
| Containers | **Docker → ghcr.io** | One image per component, tagged `build_main` |
| Data containers | **SpatialData `.zarr`** (spatial) + **AnnData `.h5ad`** (single-cell / tabular) |
| Core library | **txsim** (`theislab/txsim@dev`) | Shared spatial-transcriptomics utilities |

There are **two pipelines**: (A) dataset creation and (B) the benchmark itself.

---

## 1. Pipeline A — Dataset Creation (`process_datasets`)

Two parallel input paths (spatial + single-cell) that meet at a **combine** step. Produces the standardized inputs the benchmark consumes.

```
 SPATIAL PATH (iST)                         SINGLE-CELL PATH (scRNA-seq reference)
 ─────────────────                          ─────────────────────────────────────
 vendor raw data                            vendor raw data
   │  [spatial loader]                        │  [SC loader]  (counts only)
   ▼                                          ▼
 dataset.zarr (raw SpatialData)             dataset.h5ad (counts)
   │  [crop_region] (optional)                │  log_cp → hvg → pca → knn
   ▼                                          ▼
 cropped dataset.zarr                       processed reference.h5ad (counts + normalized + HVG/PCA/kNN)
        │                                          │
        └───────────────┬──────────────────────────┘
                        ▼
              [ process_dataset ]   ← intersect shared genes, crop huge images, rechunk
                        ▼
        raw_ist.zarr   +   scrnaseq_reference.h5ad     ← inputs to Pipeline B
```

### Loaders (`src/datasets/loaders/`)

**Spatial (iST) loaders** — raw vendor output → raw SpatialData `dataset.zarr`:
- `tenx_xenium` — 10x Xenium
- `tenx_atera` — 10x Xenium (Atera / newer format)
- `vizgen_merscope` — Vizgen MERSCOPE
- `bruker_cosmx` — Bruker CosMx
- `bruker_cosmx_nsclc` — Bruker CosMx (NSCLC)
- `allen_brain_cell_atlas_merfish` — Allen MERFISH

**Single-cell (reference) loaders** — raw → `.h5ad` (**counts only**; normalization added later):
- `allen_brain_cell_atlas`
- `andrews_human_liver_sc`
- `ganier_human_skin_sc`
- `lee_human_colon_cancer_sc`
- `lu_human_liver_cancer_sc`
- `travaglini_human_lung_sc`
- `wu_human_breast_cancer_sc`
- `zuani_human_nsclc_sc`

### Processors & shared components
- `crop_region` — optional spatial crop (local processor)
- Shared (pulled from external repos): `normalization/log_cp`, `processors/hvg | pca | knn`, `utils/extract_uns_metadata`

### Per-dataset workflows (`src/datasets/workflows/process_<name>`)
- **SC workflow:** loader → `log_cp` (adds log-CP10k `normalized` layer) → `hvg` → `pca` → `knn` → `extract_uns_metadata`. Also a generic loader-agnostic `process_scrnaseq`.
- **Spatial workflow:** loader → optional `crop_region` → `dataset.zarr`

### Combine (`process_dataset`)
Intersects spatial + SC on shared genes (`feature_name`), crops oversized images (> ~20000² px), rechunks to a uniform grid, renames `table`→`metadata`. **Does not normalize.**

---

## 2. Pipeline B — The Benchmark (`run_benchmark`)

The core preprocessing pipeline. **Data format changes partway through:** SpatialData `.zarr` through transcript assignment, then AnnData `.h5ad` from count aggregation onward.

### Core stage flow

```
 raw_ist.zarr
   │
 ① SEGMENTATION ──────────────────► segmentation.zarr
   │
 ② TRANSCRIPT ASSIGNMENT ─────────► transcript_assignments.zarr   (+ raw_ist.zarr)
   │                                        ⋮ format switch .zarr → .h5ad
 ③ COUNT AGGREGATION ─────────────► spatial_aggregated_counts.h5ad
   │
 ├─④ QC FILTER ───────────────────► spatial_qc_col.h5ad
 │
 ├─⑤ CELL VOLUME ─────────────────► cell_volumes.h5ad
   │
 ⑥ NORMALIZATION ─────────────────► spatial_normalized_counts.h5ad
   │
 ⑦ CELL TYPE ANNOTATION ──────────► spatial_with_cell_types.h5ad  (+ scrnaseq_reference.h5ad)
   │
 ⑧ EXPRESSION CORRECTION ─────────► spatial_corrected_counts.h5ad
   │
 ⑨ GENE EFFICIENCY CORRECTION ────► (corrected counts)
   │
 ⑩ DATA AGGREGATION ──────────────► spatial_processed_complete.zarr
   │
 ⑪ METRICS  ──────────────────────► score.h5ad
```

### Stage-by-stage detail

| # | Stage (namespace) | Purpose | Methods available |
|---|-------------------|---------|-------------------|
| ① | **Segmentation** (`methods_segmentation`) | Partition the image into cell label maps | `binning`, `cellpose`, `cellposev4`, `custom_segmentation`, `stardist`, `watershed` |
| ② | **Transcript assignment** (`methods_transcript_assignment`) | Assign each transcript to a cell | `basic_transcript_assignment`, `baysor`, `clustermap`, `comseg`, `fastreseg`, `pciseq`, `proseg`, `segger` |
| ③ | **Count aggregation** (`methods_count_aggregation`) | Aggregate transcripts → cell × gene counts | `basic_count_aggregation` |
| ④ | **QC filter** (`methods_qc_filter`) | Flag cells passing quality control | `basic_qc_filter` |
| ⑤ | **Cell volume** (`methods_calculate_cell_volume`) | Compute per-cell volume | `alpha_shapes` |
| ⑥ | **Normalization** (`methods_normalization`) | Normalize expression | `normalize_by_counts`, `normalize_by_volume`, `spanorm` |
| ⑦ | **Cell type annotation** (`methods_cell_type_annotation`) | Label cell types (uses scRNA-seq reference) | `mapmycells`, `moscot`, `rctd`, `singler`, `ssam`, `tacco`, `tangram` |
| ⑧ | **Expression correction** (`methods_expression_correction`) | Correct/denoise expression | `no_correction`, `resolvi_correction`, `split` |
| ⑨ | **Gene efficiency correction** (`methods_gene_efficiency_correction`) | Correct per-gene detection efficiency | `gene_efficiency_correction`, `no_correction` |
| ⑩ | **Data aggregation** (`methods_data_aggregation`) | Merge raw + processed files into one zarr | `aggregate_spatial_data` |

> Method counts: **6 segmentation · 8 transcript assignment · 7 cell-type annotation** are the "wide" stages; most others have 1–3.

---

## 3. Control Methods (`control_methods`)

Baselines that **bypass the entire pipeline** — built directly from the scRNA-seq reference, mixed in just before scoring:
- `identity` — identical copy of the scRNAseq reference (positive control / upper bound)
- `permute_celltype_annotations` — randomly permuted cell-type labels (negative control / lower bound)

Controls only produce the files the **similarity** metric needs, so **controls receive only similarity scores, not quality scores.**

---

## 4. Metrics (`metrics`)

Two metric components, each emitting several sub-metrics into `score.h5ad`.

### `quality` — reference-free (gated on data aggregation)
- `proportion_of_assigned_reads` — proportion of assigned reads
- `proportion_of_annotated_cells` — proportion of cells annotated with a cell type
- `number_of_cells` — count of cells in the spatial dataset

### `similarity` — reference-based (spatial vs. scRNA-seq reference)
- `negative_marker_purity_reads` — % negative-marker reads assigned to correct cell types
- `negative_marker_purity_cells` — % cells without negative-marker counts for their type
- `coexpr_similarity` — similarity of co-expression patterns (spatial vs. sc)
- `coexpr_similarity_celltype` — within-cell-type co-expression similarity
- `rel_pairwise_ct_expr_sim` — mean expression difference between cell-type pairs
- `rel_pairwise_gene_expr_sim` — mean expression difference between gene pairs
- `knn_mixing` — modality mixing in the joint kNN graph

---

## 5. Data Formats (nodes between stages)

| File | Container | Produced by | Contents |
|------|-----------|-------------|----------|
| `raw_ist.zarr` | SpatialData | Pipeline A | Raw imaging + transcripts (benchmark entry point) |
| `scrnaseq_reference.h5ad` | AnnData | Pipeline A | SC reference: counts + normalized + HVG/PCA/kNN |
| `segmentation.zarr` | SpatialData | ① | Cell label maps |
| `transcript_assignments.zarr` | SpatialData | ② | Per-transcript cell assignment (point cloud) |
| `spatial_aggregated_counts.h5ad` | AnnData | ③ | Raw cell × gene counts |
| `spatial_qc_col.h5ad` | AnnData | ④ | QC pass/fail column |
| `cell_volumes.h5ad` | AnnData | ⑤ | Per-cell volume column |
| `spatial_normalized_counts.h5ad` | AnnData | ⑥ | Normalized expression layer |
| `spatial_with_cell_types.h5ad` | AnnData | ⑦ | Counts + normalized + cell types |
| `spatial_corrected_counts.h5ad` | AnnData | ⑧ | Corrected + uncorrected normalized expression |
| `spatial_processed_complete.zarr` | SpatialData | ⑩ | All raw + processed data combined |
| `score.h5ad` | AnnData | ⑪ | Metric ids + values |

---

## 6. Orchestration Invariants (how the matrix runs)

Important for accurately depicting the DAG:

- **Star-matrix, NOT a full cartesian product.** At most one non-default method per run — the matrix varies **one stage at a time** around a default pipeline.
- **Controls bypass the pipeline** and are injected just before scoring.
- **Quality metric is gated** on the data-aggregation output; `aggregate_spatial_data` is its sole feeder.
- **Two normalization branches:** the standard volume path, plus an active "direct-normalization" branch (`normalize_by_counts` + `spanorm`) that skips cell-volume.
- **Direct-assignment branch is inactive** (empty / commented out).
- **Scores output:** `score_uns.yaml` (+ `dataset_uns`, `method_configs`, `metric_configs`, `task_info`, `trace.txt`) per run.

### Nextflow resource labels (for annotating compute cost)
- **Memory:** `lowmem` 20 GB · `midmem` 50 GB · `highmem` 100 GB
- **CPU:** `lowcpu` 5 · `midcpu` 15 · `highcpu` 30
- **Time:** `lowtime` 1h · `midtime` 4h · `hightime` 8h · `veryhightime` 24h
- Some methods (e.g. `segger`) are **GPU-only**.

---

## 7. Deployment Model (three gates before a run)

A component change only reaches a cluster run after **three** independent things happen — a failure at any gate means the cluster keeps using the **old container**:

```
 [1] SOURCE          →   [2] CODEGEN          →   [3] CONTAINER
 committed to            build/main deploy        Docker image pushed to
 origin/main            branch regenerated       ghcr.io @ build_main tag
                         (inlined main.nf)        (revision label = git SHA)
```

The decisive freshness signal is the image's `org.opencontainers.image.revision` label vs. `origin/main` HEAD.

---

## 8. Claude Code Skills (repo tooling, `.claude/skills/`)

Custom skills that support the *development/ops* workflow around the pipeline (not part of the data flow itself):

| Skill | What it does | When it fires |
|-------|--------------|---------------|
| **check-component** | Verifies a component is fully deployed across the 3 gates (source → codegen → container); reads the ghcr revision label without pulling the image | Before (re)running a benchmark; when an edit seems to have no effect; to confirm an image isn't stale |
| **debug-component-k8s** | Runs a component's container **live on the k8s cluster** to reproduce and fix runtime environment errors (ImportError, missing `.so`, CLI API changes, base-image issues) | When a component fails *inside* its image at runtime and needs reproduction without a full Nextflow run |
| **document-method-troubleshooting** | Captures troubleshooting into a two-layer record: in-repo `NOTES.md` (authoritative how/why/where-it-breaks) + a memory file (pointer + dated iteration log) | After debugging a component; "document this / write it up / add to the notes"; append after each further iteration |

---

## Quick element inventory (for legends / counts)

- **Pipeline stages (core):** 10 (① segmentation → ⑩ data aggregation) + metrics
- **Total method components:** ~30 across all stages
- **Widest stages:** transcript assignment (8), cell-type annotation (7), segmentation (6)
- **Control methods:** 2
- **Metric components:** 2 (quality = 3 sub-metrics, similarity = 7 sub-metrics)
- **Spatial vendors:** 4 (Xenium, MERSCOPE, CosMx, MERFISH)
- **Dataset loaders:** 6 spatial + 8 single-cell
- **Data formats:** 2 containers (`.zarr` SpatialData, `.h5ad` AnnData) across 12 file specs
- **Dev skills:** 3
