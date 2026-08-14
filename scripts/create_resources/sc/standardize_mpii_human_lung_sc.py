#!/usr/bin/env python
"""Standardize the MPII custom human-lung scRNAseq annotation into the
OpenProblems SC *loader-output* format.

Input : the pre-annotated ``current_annotation.h5ad`` (101,141 cells x 20,497
        genes, ENSEMBL var index, HLCA-style cell-type hierarchy).
Output: a clean ``dataset.h5ad`` mirroring what the SC loaders emit
        (e.g. ``wu_human_breast_cancer_sc``): a ``counts`` layer, a ``cell_type``
        hierarchy in ``.obs``, ``feature_name`` (gene symbols) in ``.var`` and
        ``dataset_*`` metadata in ``.uns``.

Deliberately NOT computed here (added by the downstream SC process step
loader -> log_cp -> hvg -> pca -> knn): ``normalized`` layer, HVG flags, PCA,
kNN graph. This file is the equivalent of the loader's ``output_raw``.

Cell-type mapping (per user choice):
    cell_type        <- celltype       (finest, 67 classes)
    cell_type_level2 <- celltype_l2    (21 classes)
    cell_type_level3 <- celltype_l3    (41 classes)
    cell_type_level4 <- celltype_l4    (67 classes)
    cell_type_level1 <- celltype_l1    (4 classes, kept for reference)

Run locally with an env that has anndata + scipy, e.g.::

    /opt/miniconda3/envs/spatialdata/bin/python \
        scripts/create_resources/sc/standardize_mpii_human_lung_sc.py \
        --input  ~/Downloads/current_annotation.h5ad \
        --output ~/Downloads/mpii_human_lung_sc/dataset.h5ad

Then upload, e.g.::

    aws s3 cp ~/Downloads/mpii_human_lung_sc/dataset.h5ad \
        s3://openproblems-data/resources/datasets/mpii_human_lung_sc/current_annotation/dataset.h5ad
"""

import argparse
import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# ── Dataset metadata ────────────────────────────────────────────────────────
DATASET_ID = "mpii_human_lung_sc"
DATASET_NAME = "MPII human lung scRNAseq reference (current annotation)"
DATASET_URL = ""  # custom / unpublished MPII reference
DATASET_REFERENCE = ""
DATASET_SUMMARY = (
    "Custom MPII human lung scRNAseq reference with a fine-grained "
    "(67 cell types) HLCA-style annotation hierarchy."
)
DATASET_DESCRIPTION = (
    "Custom single-cell RNA-seq reference of healthy human lung generated at MPII, "
    "used as the annotation reference for the MPII Xenium human lung dataset. "
    "101,141 cells x 20,497 genes with a 4-level cell-type hierarchy (celltype_l1..l4)."
)
DATASET_ORGANISM = "homo_sapiens"

# standardized cell-type column  <-  source column in current_annotation.h5ad
CELLTYPE_MAP = {
    "cell_type": "celltype",        # primary label (finest, == celltype_l4) — user choice
    "cell_type_level1": "celltype_l1",
    "cell_type_level2": "celltype_l2",
    "cell_type_level3": "celltype_l3",
    "cell_type_level4": "celltype_l4",
}

# extra metadata carried over from .obs if present (standardized -> source)
EXTRA_OBS = {
    "donor_id": "mpii_subject_id",
    "sample": "sample",
    "assay": "assay",
    "sex": "sex",
    "disease": "disease",
    "tissue": "tissue",
    "organism": "organism",
}

_MISSING = {"nan", "None", "NA", "NaN", "", "unknown", "Unknown"}


def _symbols_from_var(var: pd.DataFrame, index: pd.Index) -> np.ndarray:
    """Gene symbols for .var: prefer gene_names/gene_name, fall back to ENSEMBL id."""
    symbols = None
    for col in ("gene_names", "gene_name"):
        if col in var.columns:
            symbols = pd.Series(var[col].astype(str).values, index=range(len(var)))
            break
    ensembl = (
        var["gene_ids"].astype(str).values
        if "gene_ids" in var.columns
        else index.astype(str).values
    )
    if symbols is None:
        symbols = pd.Series(index.astype(str).values, index=range(len(var)))
    # never leave a null/placeholder feature_name — fall back to the ENSEMBL id
    symbols = symbols.replace({m: np.nan for m in {"nan", "None", "NA", "NaN", ""}})
    symbols = symbols.fillna(pd.Series(ensembl, index=range(len(var))))
    return symbols.values, ensembl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input",
        default=str(Path.home() / "Downloads" / "current_annotation.h5ad"),
        help="Path to current_annotation.h5ad",
    )
    ap.add_argument(
        "--output",
        default=str(Path.home() / "Downloads" / "mpii_human_lung_sc" / "dataset.h5ad"),
        help="Path to write the standardized dataset.h5ad",
    )
    ap.add_argument(
        "--counts-from",
        choices=["X", "raw"],
        default="X",
        help="Which matrix holds the raw counts (both are integer in this file; X is aligned to .var/.obs).",
    )
    args = ap.parse_args()

    print(f"Reading {args.input} ...", flush=True)
    adata = ad.read_h5ad(args.input)
    print(f"  loaded: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    # ── counts ──────────────────────────────────────────────────────────────
    if args.counts_from == "raw":
        if adata.raw is None:
            raise SystemExit("--counts-from raw requested but .raw is absent")
        counts = adata.raw.X
        var_src = adata.raw.var
    else:
        counts = adata.X
        var_src = adata.var
    counts = counts.tocsr() if hasattr(counts, "tocsr") else csr_matrix(counts)
    counts = counts.astype("float32")
    if counts.min() < 0:
        raise SystemExit("counts contain negative values — X is likely not raw counts; try --counts-from raw")

    # ── var: feature_name (symbols) + feature_id (ENSEMBL) ────────────────────
    symbols, ensembl = _symbols_from_var(var_src, adata.var_names)
    var = pd.DataFrame(index=adata.var_names.astype(str))
    var["feature_name"] = symbols
    var["gene_symbol"] = symbols
    var["feature_id"] = ensembl

    # ── obs: cell_type hierarchy + carried-over metadata ─────────────────────
    obs = pd.DataFrame(index=adata.obs_names.astype(str))
    for new, src in CELLTYPE_MAP.items():
        if src in adata.obs.columns:
            obs[new] = adata.obs[src].astype(str).values
        else:
            print(f"  WARNING: source obs column '{src}' not found; skipping '{new}'")
    for new, src in EXTRA_OBS.items():
        if src in adata.obs.columns and new not in obs.columns:
            obs[new] = adata.obs[src].astype(str).values

    if "cell_type" not in obs.columns:
        raise SystemExit("primary source column 'celltype' not found — cannot set cell_type")

    # free the large original (incl. .raw, obsm, obsp) before building the output
    del adata
    gc.collect()

    # ── assemble clean AnnData ───────────────────────────────────────────────
    new = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    new.layers["counts"] = counts

    # similarity metric forbids NaN / "None" in cell_type — drop those cells
    bad = new.obs["cell_type"].str.strip().str.lower().isin({m.lower() for m in _MISSING})
    bad = bad | new.obs["cell_type"].isna()
    if bool(bad.any()):
        print(f"  dropping {int(bad.sum())} cells with missing/placeholder cell_type", flush=True)
        new = new[~bad.values].copy()

    # cell-type columns must be categorical (similarity reads .dtype.categories)
    for col in [c for c in new.obs.columns if c.startswith("cell_type")]:
        new.obs[col] = pd.Categorical(new.obs[col].astype(str))

    # ── dataset metadata ─────────────────────────────────────────────────────
    new.uns["dataset_id"] = DATASET_ID
    new.uns["dataset_name"] = DATASET_NAME
    new.uns["dataset_url"] = DATASET_URL
    new.uns["dataset_reference"] = DATASET_REFERENCE
    new.uns["dataset_summary"] = DATASET_SUMMARY
    new.uns["dataset_description"] = DATASET_DESCRIPTION
    new.uns["dataset_organism"] = DATASET_ORGANISM

    out = Path(args.output).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"Writing {out}\n"
        f"  {new.n_obs} cells x {new.n_vars} genes | layers={list(new.layers)} | "
        f"cell_type classes={new.obs['cell_type'].cat.categories.size}",
        flush=True,
    )
    new.write_h5ad(out, compression="gzip")
    print("Done.", flush=True)
    print(
        "\nNext steps (this file is the RAW, counts-only reference):\n"
        f"  1. Upload:  aws s3 cp {out} "
        "s3://openproblems-data/resources/raw_data/txSim_custom/MPII/mpii_human_lung_sc_raw.h5ad\n"
        "  2. Process: run scripts/create_resources/sc/process_mpii_human_lung_sc_nebius.sh\n"
        "     (generic datasets/workflows/process_scrnaseq: log_cp -> hvg -> pca -> knn).\n"
        "     It adds the 'normalized' layer, HVG, PCA and kNN and publishes the final\n"
        f"     input_sc to s3://openproblems-data/resources/datasets/{DATASET_ID}/current_annotation/dataset.h5ad\n"
        "  3. Combine: run scripts/create_resources/combine/process_datasets_mpii_nebius.sh"
    )


if __name__ == "__main__":
    main()
