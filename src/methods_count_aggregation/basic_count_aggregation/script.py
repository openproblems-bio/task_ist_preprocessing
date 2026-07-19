import numpy as np
import pandas as pd
import anndata as ad
import spatialdata as sd
from scipy.sparse import csr_matrix


def generate_adata(input_spots, cell_id_col="cell", gene_col="Gene"):
    """Aggregate per-transcript assignments into a cell x gene AnnData.

    Reimplementation of txsim.preprocessing.generate_adata: txsim's version is
    incompatible with anndata>=0.12 (it passes the removed `dtype=` argument to
    AnnData() and uses per-row item assignment `adata[cell, :] = ...`, which
    anndata no longer supports). This fills the count matrix directly and
    produces an identical output structure (X, layers['raw_counts'],
    obs/var stats, uns['spots'/'pct_noise']).
    """
    spots = input_spots.copy()
    pct_noise = sum(spots[cell_id_col] <= 0) / len(spots[cell_id_col])
    spots_raw = spots.copy()  # kept in uns['spots']; 0 (background) -> None
    spots_raw.loc[spots_raw[cell_id_col] == 0, cell_id_col] = None
    spots = spots[spots[cell_id_col] > 0]

    cell_ids = pd.unique(spots[cell_id_col])
    genes = pd.unique(spots[gene_col])
    cell_pos = {c: i for i, c in enumerate(cell_ids)}

    # Populate the count matrix + centroids per cell (no AnnData item assignment).
    # feature_name may be categorical, so value_counts() can return unobserved
    # categories; reindex to `genes` (fill 0) to align to the var order, mirroring
    # txsim's `.reindex(var_names, fill_value=0)`.
    X = np.zeros((len(cell_ids), len(genes)), dtype=np.float32)
    centroid_x = np.zeros(len(cell_ids))
    centroid_y = np.zeros(len(cell_ids))
    for cell_id, grp in spots.groupby(cell_id_col, sort=False):
        row = cell_pos[cell_id]
        cts = grp[gene_col].value_counts().reindex(genes, fill_value=0)
        X[row, :] = cts.values
        centroid_x[row] = grp["x"].mean()
        centroid_y[row] = grp["y"].mean()

    adata = ad.AnnData(csr_matrix(X))
    adata.obs["cell_id"] = cell_ids
    adata.obs_names = [f"{i:d}" for i in cell_ids]
    adata.var_names = genes
    adata.obs["centroid_x"] = centroid_x
    adata.obs["centroid_y"] = centroid_y

    adata.uns["spots"] = spots_raw
    adata.uns["pct_noise"] = pct_noise
    adata.layers["raw_counts"] = adata.X.copy()

    adata.obs["n_counts"] = np.ravel(adata.layers["raw_counts"].sum(axis=1))
    adata.obs["n_genes"] = adata.layers["raw_counts"].getnnz(axis=1)
    adata.var["n_counts"] = np.ravel(adata.layers["raw_counts"].sum(axis=0))
    adata.var["n_cells"] = adata.layers["raw_counts"].getnnz(axis=0)
    return adata


## VIASH START
par = {
  'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr',
  'output': 'raw_counts.h5ad',
}
meta = {
  'name': 'basic'
}
## VIASH END

sdata = sd.read_zarr(par['input'])
df = sdata['transcripts'].compute() # TODO: Could optimize tx.preprocessing.generate_adata to work on spatialdata

adata = generate_adata(df, cell_id_col='cell_id', gene_col='feature_name') #TODO: x and y refers to a specific coordinate system. Decide which space we want to use here. (probably should be handled in the previous assignment step)
adata.layers['counts'] = adata.layers['raw_counts']
del adata.layers['raw_counts']
adata.var["gene_name"] = adata.var_names

# currently the function also saves the transcripts in the adata object, but this is not necessary here
del adata.uns['spots']
del adata.uns['pct_noise']

adata.write_h5ad(par['output'], compression="gzip")
