import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.sparse import issparse
from scanpy._utils import view_to_actual
from sklearn.utils import sparsefuncs

## VIASH START
par = {
  'input_spatial_aggregated_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad',
  'input_cell_volumes': 'resources_test/task_ist_preprocessing/mouse_brain_combined/cell_volumes.h5ad',
  'output': 'spatial_norm_counts.h5ad',

  'target_volume': None,
}
meta = {
  'name': 'normalize_by_volume',
}
## VIASH END

assert par['input_cell_volumes'] is not None, 'Volume input is required for this normalization method.'


### Function definition ###

def normalize_by_volume(
    adata: ad.AnnData,
    target_volume: float | None = None,
    volume_key: str = 'volume',
    key_added: str | None = None,
    layer: str | None = None,
    inplace: bool = True,
) -> np.ndarray | sparse.spmatrix | None:
    """Normalize counts by cell volume
    
    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    target_volume
        If `None`, after normalization, each observation (cell) has a
        volume equal to the median of cell volumes. Note that when providing
        a float value it is crucial to know in which spatial units the cell volume 
        was calculated.
    volume_key : str
        Name of the field in `adata.obs` where the cell volume (or area) is stored, 
        by default 'volume'
    key_added
        Name of the field in `adata.obs` where the normalization factor is
        stored.
    layer
        Layer to normalize instead of `X`. If `None`, `X` is normalized.
    inplace
        Whether to update `adata` or return dictionary with normalized copies of
        `adata.X` and `adata.layers`.
    
    
    Returns
    -------
    np.ndarray
        If ``inplace=True``, ``adata.X`` is updated with the normalized values. 
        Otherwise, returns the normalized numpy array or sparse matrix
    """
    
    view_to_actual(adata)
    
    X = adata.X if (layer is None) else adata.layers[layer]
    
    assert volume_key in adata.obs.columns, f"key {volume_key} not found in adata.obs"
    assert not adata.obs[volume_key].isnull().any(), "Nan values found in adata.obs[volume_key]. All cells must have valid volumes."
    assert (adata.obs[volume_key] > 0).all(), "All cell volumes in adata.obs[volume_key] must be > 0"
    
    if target_volume is None:
        target_volume = adata.obs[volume_key].median()
    
    volumes = adata.obs[volume_key].values.copy()
    norm_factors = target_volume / volumes
    
    if issparse(X):
        if X.dtype != norm_factors.dtype:
            X = X.astype(np.float64)
            norm_factors = norm_factors.astype(np.float64)
        if inplace:
            sparsefuncs.inplace_row_scale(X, norm_factors)
        else:
            X_norm = X.copy()
            sparsefuncs.inplace_row_scale(X_norm, norm_factors)
    else:
        if inplace:
            np.divide(X, volumes[:,None], out=X)
        else:
            X_norm = X / volumes[:, None]
    
    if key_added:
        adata.obs[key_added] = norm_factors
    
    if not inplace:
        return X_norm


### Main function ###

print('Reading input files', flush=True)
adata = ad.read_h5ad(par['input_spatial_aggregated_counts'])
adata_volume = ad.read_h5ad(par['input_cell_volumes'])
assert adata.obs["cell_id"].astype(int).sort_values().equals(adata_volume.obs["cell_id"].astype(int).sort_values()), "Cell IDs do not match"
adata.obs["volume"] = adata_volume.obs.loc[adata.obs["cell_id"].astype(str), "volume"]

print('Normalizing by volume', flush=True)
adata.layers['normalized'] = adata.layers['counts'].copy()
normalize_by_volume(adata, layer="normalized", target_volume=par["target_volume"])
sc.pp.log1p(adata, layer="normalized")

print('Writing output', flush=True)
adata.write(par['output'])