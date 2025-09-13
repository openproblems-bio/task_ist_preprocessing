import numpy as np
import pandas as pd
import anndata as ad
import spatialdata as sd

def create_dummy_transcript_assignment_table(adata: ad.AnnData) -> sd.SpatialData:
    """ Create a dummy transcript assignment table from an AnnData object.

    Arguments
    ---------
    adata: ad.AnnData
        The AnnData object to create a dummy transcript assignment table from.

    Returns
    -------
    sdata_transcripts_only: sd.SpatialData
        The SpatialData object with the dummy transcript assignment table.
    """

    # Convert the sparse matrix to coo for access of row and col as arrays
    coo = adata.layers["counts"].tocoo()

    # Get cell and gene vectors
    counts = np.astype(coo.data, np.int64)
    cell_id_idx = np.repeat(coo.row, counts)
    gene_id_idx = np.repeat(coo.col, counts)

    obs_names = adata.obs_names.values
    var_names = adata.var_names.values

    cell_ids = obs_names[cell_id_idx]
    genes = var_names[gene_id_idx]

    # Create the dummy transcript assignment table
    transcripts_df = pd.DataFrame({
        "x": cell_id_idx,
        "y": cell_id_idx,
        "z": 0,
        "feature_name": genes,
        "cell_id": cell_ids,
        "overlaps_nucleus": 0,
        "qv": 0,
        "transcript_id": [i for i in range(len(cell_ids))]
    })

    # Create the transcripts sdata
    sdata_table = ad.AnnData(obs=adata.obs[[]], var=adata.var[[]])
    sdata_table.obs["cell_id"] = adata.obs_names.values

    sdata_transcripts_only = sd.SpatialData(
        points={"transcripts": sd.models.PointsModel.parse(transcripts_df)},
        tables={"table": sdata_table}
    )
    
    return sdata_transcripts_only


def add_layers_obs_var_to_scrnaseq_ref(adata: ad.AnnData) -> ad.AnnData:
    """ Add layers, obs and var columns to an AnnData object to have the same structure as the processed spatial data.

    Arguments
    ---------
    adata: ad.AnnData
        The AnnData object to add layers, obs and var columns to.
    """

    adata.layers["normalized_uncorrected"] = adata.layers["normalized"]
    adata.obs["cell_id"] = adata.obs.index
    adata.obs["centroid_x"] = 0
    adata.obs["centroid_y"] = 0
    adata.obs["centroid_z"] = 0
    adata.obs["n_counts"] = np.array(adata.layers["counts"].sum(axis=1))[:,0]
    adata.obs["n_genes"] = np.array((adata.layers["counts"] > 0).sum(axis=1))[:,0]
    adata.obs["volume"] = 1
    adata.var["gene_name"] = adata.var.index
    adata.var["n_counts"] = np.array(adata.layers["counts"].sum(axis=0))[0,:]
    adata.var["n_cells"] = np.array((adata.layers["counts"] > 0).sum(axis=0))[0,:]