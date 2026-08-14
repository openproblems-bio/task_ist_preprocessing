#!/usr/bin/env python
"""Standardize the custom LTX human-lung scRNAseq annotation into the
OpenProblems SC *loader-output* format (raw, counts-only).

Input : ltx_full_annotation_2026_04_25.h5ad (363,099 cells x 18,117 genes;
        gene SYMBOLS in var index, ENSEMBL in var['gene_ids']; raw counts in
        layers['counts']; 25-class annotation in obs['celltype_final']).
Output: a clean dataset.h5ad mirroring what the SC loaders emit (a `counts`
        layer, `cell_type` in .obs, `feature_name` (symbols) in .var, and
        `dataset_*` metadata in .uns).

Deliberately NOT computed here (added downstream by process_scrnaseq:
log_cp -> hvg -> pca -> knn): the `normalized` layer, HVG flags, PCA, kNN.

Memory: the source file also carries X, a normcounts layer and a large obsp
(neighbour graph). We therefore read ONLY layers['counts'] + obs + var via
anndata.io.read_elem, so peak memory ~= the counts matrix, not the full 20 GB.

Run locally with an env that has anndata>=0.11 + scipy, e.g.::

    /opt/miniconda3/envs/spatialdata/bin/python \
        scripts/create_resources/sc/standardize_ltx_sc.py \
        --input  ~/projects/tests_fot_txsim/ltx_full_annotation_2026_04_25.h5ad \
        --output ~/projects/tests_fot_txsim/ltx_human_lung_sc/dataset.h5ad

Then upload (profile 'op')::

    aws s3 cp --profile op <output> \
        s3://openproblems-data/resources/raw_data/txSim_custom/ltx/ltx_human_lung_sc_raw.h5ad
"""

import argparse
import gc
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

# ── Dataset metadata ────────────────────────────────────────────────────────
DATASET_ID = "ltx_human_lung_sc"
DATASET_NAME = "LTX human lung scRNAseq reference (full annotation 2026-04-25)"
DATASET_URL = ""            # custom / unpublished LTX reference
DATASET_REFERENCE = ""
DATASET_SUMMARY = (
    "Custom LTX human lung scRNAseq reference with a 25-class cell-type "
    "annotation (celltype_final), used as the annotation reference for the "
    "LTX Xenium human lung samples."
)
DATASET_DESCRIPTION = (
    "Single-cell RNA-seq reference of human lung (LTX cohort), used as the "
    "annotation reference for the LTX Xenium dataset. 363,099 cells x 18,117 "
    "genes; primary label in celltype_final (25 classes)."
)
DATASET_ORGANISM = "homo_sapiens"

# primary label (per user): cell_type <- celltype_final. The coarser LTX columns
# ('lineage', 'celltype coarse', 'celltype refined') could optionally be mapped
# to cell_type_level2..4 (all optional in file_common_scrnaseq); left out here to
# avoid guessing granularity — add if the benchmark needs a hierarchy.
CELLTYPE_COL = "celltype_final"
COUNTS_LAYER = "counts"

# optional metadata carried over from .obs if present (standardized -> source)
EXTRA_OBS = {
    "donor_id": "sampleID",
    "sample": "Sample Pseudonym",
    "sample_type": "sample type",
    "batch": "batch",
    "lineage": "lineage",
}

_MISSING = {"nan", "none", "na", "", "unknown"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Path to the raw LTX annotation h5ad")
    ap.add_argument("--output", required=True, help="Path to write the standardized dataset.h5ad")
    args = ap.parse_args()

    inp = Path(args.input).expanduser()
    print(f"Reading (selectively) {inp} ...", flush=True)
    with h5py.File(inp, "r") as f:
        if COUNTS_LAYER not in f.get("layers", {}):
            raise SystemExit(f"layers['{COUNTS_LAYER}'] not found in {inp}")
        counts = read_elem(f["layers"][COUNTS_LAYER])
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
    print(f"  loaded counts {counts.shape} + obs{obs.shape} + var{var.shape}", flush=True)

    # ── counts: enforce raw, non-negative, csr float32 ───────────────────────
    counts = counts.tocsr() if hasattr(counts, "tocsr") else csr_matrix(counts)
    counts = counts.astype("float32")
    if counts.min() < 0:
        raise SystemExit(f"layers['{COUNTS_LAYER}'] contains negative values — not raw counts")
    frac_int = float(np.mean(counts.data == np.round(counts.data))) if counts.nnz else 1.0
    print(f"  counts: nnz={counts.nnz:,} min={counts.min()} max={counts.max()} integer-fraction={frac_int:.3f}", flush=True)
    # Enforce raw integer counts. Downstream count-based methods (RCTD, SPLIT) call
    # spacexr::check_counts, which hard-rejects non-integer references. A processed
    # atlas export can carry log-normalized values in layers['counts'] (and in .raw.X)
    # while retaining only per-cell totals in obs — do NOT let that pass silently as a
    # "reference". If this fires, point COUNTS_LAYER at the true raw count matrix.
    if frac_int < 0.999:
        raise SystemExit(
            f"layers['{COUNTS_LAYER}'] is not raw integer counts "
            f"(integer-fraction={frac_int:.3f}, max={counts.max():.4g}) — looks normalized/"
            f"log-transformed. RCTD/SPLIT require integer counts; standardize needs the raw "
            f"count matrix, not this layer."
        )

    # ── var: feature_name (symbols, = index) + feature_id (ENSEMBL) ──────────
    var_index = var.index.astype(str)
    new_var = pd.DataFrame(index=var_index)
    new_var["feature_name"] = var_index.values
    new_var["gene_symbol"] = var_index.values
    new_var["feature_id"] = (
        var["gene_ids"].astype(str).values if "gene_ids" in var.columns else var_index.values
    )

    # ── obs: cell_type (from celltype_final) + carried metadata ──────────────
    if CELLTYPE_COL not in obs.columns:
        raise SystemExit(f"obs['{CELLTYPE_COL}'] not found — cannot set cell_type")
    new_obs = pd.DataFrame(index=obs.index.astype(str))
    new_obs["cell_type"] = obs[CELLTYPE_COL].astype(str).values
    for new, src in EXTRA_OBS.items():
        if src in obs.columns:
            new_obs[new] = obs[src].astype(str).values

    del obs, var
    gc.collect()

    # ── assemble clean AnnData (X and layers['counts'] share the matrix) ─────
    adata = ad.AnnData(X=counts, obs=new_obs, var=new_var)
    adata.layers["counts"] = counts

    # similarity metric forbids NaN / "None" in cell_type — drop those cells
    bad = adata.obs["cell_type"].str.strip().str.lower().isin(_MISSING) | adata.obs["cell_type"].isna()
    if bool(bad.any()):
        print(f"  dropping {int(bad.sum())} cells with missing/placeholder cell_type", flush=True)
        adata = adata[~bad.values].copy()

    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"].astype(str))

    # ── dataset metadata ─────────────────────────────────────────────────────
    adata.uns["dataset_id"] = DATASET_ID
    adata.uns["dataset_name"] = DATASET_NAME
    adata.uns["dataset_url"] = DATASET_URL
    adata.uns["dataset_reference"] = DATASET_REFERENCE
    adata.uns["dataset_summary"] = DATASET_SUMMARY
    adata.uns["dataset_description"] = DATASET_DESCRIPTION
    adata.uns["dataset_organism"] = DATASET_ORGANISM

    out = Path(args.output).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"Writing {out}\n"
        f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes | layers={list(adata.layers)} | "
        f"cell_type classes={adata.obs['cell_type'].cat.categories.size}",
        flush=True,
    )
    adata.write_h5ad(out, compression="gzip")
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
