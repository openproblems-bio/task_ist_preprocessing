#!/usr/bin/env python
"""Standardize a Kuppe kidney Flex-snRNA ``annotated.h5ad`` into the OpenProblems
SC *loader-output* format (raw, counts-only).

The Kuppe kidney MERFISH release ships two snRNA-seq references under
``Flex-snRNA/{Cancer,Healthy}/annotated.h5ad``. They use DIFFERENT cell-type
vocabularies and gene sets, so this script standardizes ONE file per run — point
it at each in turn to get two independent references (as chosen for this dataset).

Input  (per file): ``annotated.h5ad`` with
        - raw integer counts in ``layers['counts']`` (``X`` / ``layers['log1p_norm']``
          are already normalized — ignored here; the ``normalized`` layer is
          recomputed downstream),
        - gene SYMBOLS in ``var.index`` and ENSEMBL ids in ``var['gene_ids']``,
        - the curated cell-type label in ``obs['annotation']`` (Cancer: 13 classes
          incl. ``Tumor`` and some empty-string entries; Healthy: 19 kidney classes).
          ``celltype_hint`` / ``predictions_unconstrained`` are automated ontology
          hints, NOT the curated label — deliberately not used.
Output: a clean ``dataset.h5ad`` mirroring what the SC loaders emit (a ``counts``
        layer, ``cell_type`` in ``.obs``, ``feature_name`` (symbols) in ``.var`` and
        ``dataset_*`` metadata in ``.uns``).

Deliberately NOT computed here (added downstream by process_scrnaseq:
log_cp -> hvg -> pca -> knn): the ``normalized`` layer, HVG flags, PCA, kNN.

Memory: the source files also carry ``X``, a ``log1p_norm`` layer and large
``obsp``/``uns`` (neighbour graph, leiden markers). We therefore read ONLY
``layers['counts']`` + ``obs`` + ``var`` via ``anndata.io.read_elem``, so peak
memory ~= the counts matrix.

Run locally with an env that has anndata>=0.11 + scipy, e.g. (both conditions)::

    PY=/opt/miniconda3/envs/spatialdata/bin/python
    BASE="/Volumes/SeagateHHD/Kuppe_kidney_merfish/240824_download/Flex-snRNA"

    $PY scripts/create_resources/sc/standardize_kuppe_kidney_sc.py \
        --condition cancer \
        --input  "$BASE/Cancer/annotated.h5ad" \
        --output ~/Downloads/kuppe_kidney_cancer_sc/dataset.h5ad

    $PY scripts/create_resources/sc/standardize_kuppe_kidney_sc.py \
        --condition healthy \
        --input  "$BASE/Healthy/annotated.h5ad" \
        --output ~/Downloads/kuppe_kidney_healthy_sc/dataset.h5ad

Then upload (profile 'op'), matching the paths in process_kuppe_kidney_sc_nebius.sh::

    aws s3 cp --profile op ~/Downloads/kuppe_kidney_cancer_sc/dataset.h5ad \
        s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney/cancer/dataset.h5ad
    aws s3 cp --profile op ~/Downloads/kuppe_kidney_healthy_sc/dataset.h5ad \
        s3://openproblems-data/resources/raw_data/txSim_custom/kuppe_kidney/healthy/dataset.h5ad
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

# ── Per-condition metadata presets ──────────────────────────────────────────
# One reference per condition (chosen because the two files use disjoint cell-type
# vocabularies). Every field is overridable on the CLI.
PRESETS = {
    "cancer": dict(
        dataset_id="kuppe_kidney_cancer_sc",
        dataset_name="Kuppe kidney cancer snRNA-seq reference (Flex)",
        dataset_summary=(
            "Kuppe lab kidney cancer single-nucleus RNA-seq reference (10x Flex), "
            "used as the annotation reference for the Kuppe kidney MERFISH dataset."
        ),
        dataset_description=(
            "Single-nucleus RNA-seq reference of kidney cancer tissue (Kuppe lab, 10x "
            "Flex), curated cell-type label in obs['annotation'] (incl. a Tumor class). "
            "Used to annotate the Kuppe kidney MERFISH spatial data."
        ),
        disease="kidney cancer",
    ),
    "healthy": dict(
        dataset_id="kuppe_kidney_healthy_sc",
        dataset_name="Kuppe kidney healthy snRNA-seq reference (Flex)",
        dataset_summary=(
            "Kuppe lab healthy kidney single-nucleus RNA-seq reference (10x Flex), "
            "used as the annotation reference for the Kuppe kidney MERFISH dataset."
        ),
        dataset_description=(
            "Single-nucleus RNA-seq reference of healthy kidney tissue (Kuppe lab, 10x "
            "Flex), curated kidney cell-type label in obs['annotation']. Used to annotate "
            "the Kuppe kidney MERFISH spatial data."
        ),
        disease="normal",
    ),
}

DATASET_ORGANISM = "homo_sapiens"

# placeholder/missing cell-type strings to drop (the similarity metric forbids
# NaN/"None" in cell_type; the Cancer file has empty-string annotations)
_MISSING = {"nan", "none", "na", "", "unknown"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Path to the Kuppe kidney annotated.h5ad")
    ap.add_argument("--output", required=True, help="Path to write the standardized dataset.h5ad")
    ap.add_argument(
        "--condition",
        choices=sorted(PRESETS),
        help="Fills dataset metadata + disease from a preset (cancer/healthy). Optional if all --dataset-* are given.",
    )
    ap.add_argument("--celltype-col", default="annotation", help="obs column holding the curated cell type (default: annotation)")
    ap.add_argument("--counts-layer", default="counts", help="layer holding raw integer counts (default: counts)")
    # metadata overrides (default to the preset when --condition is given)
    ap.add_argument("--dataset-id")
    ap.add_argument("--dataset-name")
    ap.add_argument("--dataset-summary")
    ap.add_argument("--dataset-description")
    ap.add_argument("--dataset-url", default="")
    ap.add_argument("--dataset-reference", default="")
    ap.add_argument("--disease")
    ap.add_argument("--organism", default=DATASET_ORGANISM)
    args = ap.parse_args()

    preset = PRESETS.get(args.condition, {})

    # CLI override wins, else preset value.
    dataset_id = args.dataset_id or preset.get("dataset_id")
    dataset_name = args.dataset_name or preset.get("dataset_name")
    dataset_summary = args.dataset_summary or preset.get("dataset_summary")
    dataset_description = args.dataset_description or preset.get("dataset_description")
    disease = args.disease or preset.get("disease")
    if not (dataset_id and dataset_name and dataset_summary and dataset_description):
        raise SystemExit(
            "Missing dataset metadata: pass --condition {cancer,healthy} or all of "
            "--dataset-id/--dataset-name/--dataset-summary/--dataset-description."
        )

    inp = Path(args.input).expanduser()
    print(f"Reading (selectively) {inp} ...", flush=True)
    with h5py.File(inp, "r") as f:
        layers = f.get("layers", {})
        if args.counts_layer not in layers:
            raise SystemExit(f"layers['{args.counts_layer}'] not found in {inp} (have: {list(layers.keys())})")
        counts = read_elem(f["layers"][args.counts_layer])
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
    print(f"  loaded counts {counts.shape} + obs{obs.shape} + var{var.shape}", flush=True)

    # ── counts: enforce raw, non-negative, integer, csr float32 ──────────────
    counts = counts.tocsr() if hasattr(counts, "tocsr") else csr_matrix(counts)
    counts = counts.astype("float32")
    if counts.min() < 0:
        raise SystemExit(f"layers['{args.counts_layer}'] contains negative values — not raw counts")
    frac_int = float(np.mean(counts.data == np.round(counts.data))) if counts.nnz else 1.0
    print(f"  counts: nnz={counts.nnz:,} min={counts.min()} max={counts.max()} integer-fraction={frac_int:.3f}", flush=True)
    # Downstream count-based methods (RCTD, SPLIT) call spacexr::check_counts, which
    # hard-rejects non-integer references. The Kuppe X/log1p_norm layers ARE normalized;
    # only layers['counts'] is raw — guard against pointing at the wrong layer.
    if frac_int < 0.999:
        raise SystemExit(
            f"layers['{args.counts_layer}'] is not raw integer counts "
            f"(integer-fraction={frac_int:.3f}, max={counts.max():.4g}) — looks normalized. "
            f"Point --counts-layer at the raw count matrix."
        )

    # ── var: feature_name (symbols, = index) + feature_id (ENSEMBL) ──────────
    var_index = var.index.astype(str)
    new_var = pd.DataFrame(index=var_index)
    new_var["feature_name"] = var_index.values
    new_var["gene_symbol"] = var_index.values
    new_var["feature_id"] = (
        var["gene_ids"].astype(str).values if "gene_ids" in var.columns else var_index.values
    )

    # ── obs: cell_type (from the curated annotation) + disease ───────────────
    if args.celltype_col not in obs.columns:
        raise SystemExit(f"obs['{args.celltype_col}'] not found — cannot set cell_type (have: {list(obs.columns)[:30]})")
    new_obs = pd.DataFrame(index=obs.index.astype(str))
    new_obs["cell_type"] = obs[args.celltype_col].astype(str).values
    if disease:
        new_obs["disease"] = disease

    del obs, var
    gc.collect()

    # ── assemble clean AnnData (X and layers['counts'] share the matrix) ─────
    adata = ad.AnnData(X=counts, obs=new_obs, var=new_var)
    adata.layers["counts"] = counts

    # similarity metric forbids NaN / "None" / empty in cell_type — drop those cells
    bad = adata.obs["cell_type"].str.strip().str.lower().isin(_MISSING) | adata.obs["cell_type"].isna()
    if bool(bad.any()):
        print(f"  dropping {int(bad.sum())} cells with missing/placeholder cell_type", flush=True)
        adata = adata[~bad.values].copy()

    # cell_type must be categorical (similarity reads .dtype.categories)
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"].astype(str))
    if "disease" in adata.obs.columns:
        adata.obs["disease"] = pd.Categorical(adata.obs["disease"].astype(str))

    # ── dataset metadata ─────────────────────────────────────────────────────
    adata.uns["dataset_id"] = dataset_id
    adata.uns["dataset_name"] = dataset_name
    adata.uns["dataset_url"] = args.dataset_url
    adata.uns["dataset_reference"] = args.dataset_reference
    adata.uns["dataset_summary"] = dataset_summary
    adata.uns["dataset_description"] = dataset_description
    adata.uns["dataset_organism"] = args.organism

    out = Path(args.output).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"Writing {out}\n"
        f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes | layers={list(adata.layers)} | "
        f"cell_type classes={adata.obs['cell_type'].cat.categories.size} | dataset_id={dataset_id}",
        flush=True,
    )
    adata.write_h5ad(out, compression="gzip")
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
