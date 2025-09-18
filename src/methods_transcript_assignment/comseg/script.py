import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

## VIASH START
par = {
    "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
    "output": "transcripts.zarr",

    "transcripts_key": "transcripts",
    "shapes_key": "cell_boundaries",
    "images_key": "morphology_mip",
    "patch_width": 1200,
    "patch_overlap": 50,
    "transcript_patch_width": 200,
    "mean_cell_diameter": 15.0,
    "max_cell_radius": 25.0,
    "alpha": 0.5,
    "min_rna_per_cell": 5,
    "gene_column": "feature_name",
    "norm_vector": False,
    "allow_disconnected_polygon": True,
}
## VIASH END

def fixed_count_transcripts_aligned(geo_df, points, value_key):
    """
    The same function as sopa.aggregation.transcripts._count_transcripts_aligned.
    Minor change just the matrix X is converted to csr_matrix, to avoid bug error in comseg call

    """
    from scipy.sparse import csr_matrix
    from anndata import AnnData
    from dask.diagnostics import ProgressBar
    from functools import partial
    from sopa._settings import settings
    import geopandas as gpd
    def _add_csr(X_partitions, geo_df, partition, gene_column, gene_names ):
        if settings.gene_exclude_pattern is not None:
            partition = partition[~partition[gene_column].str.match(settings.gene_exclude_pattern, case=False, na=False)]

        points_gdf = gpd.GeoDataFrame(partition, geometry=gpd.points_from_xy(partition["x"], partition["y"]))
        joined = geo_df.sjoin(points_gdf)
        cells_indices, column_indices = joined.index, joined[gene_column].cat.codes
        cells_indices = cells_indices[column_indices >= 0]
        column_indices = column_indices[column_indices >= 0]
        X_partition = csr_matrix((np.full(len(cells_indices), 1), (cells_indices, column_indices)),
            shape=(len(geo_df), len(gene_names)),
        )
        X_partitions.append(X_partition)
    

    points[value_key] = points[value_key].astype("category").cat.as_known()
    gene_names = points[value_key].cat.categories.astype(str)
    X = csr_matrix((len(geo_df), len(gene_names)), dtype=int)
    adata = AnnData(X=X, var=pd.DataFrame(index=gene_names))
    adata.obs_names = geo_df.index.astype(str)
    geo_df = geo_df.reset_index()
    X_partitions = []
    with ProgressBar():
        points.map_partitions(
            partial(_add_csr, X_partitions, geo_df, gene_column=value_key, gene_names=gene_names),
            meta=(),
        ).compute()
    for X_partition in X_partitions:
        adata.X += X_partition
    if settings.gene_exclude_pattern is not None:
        adata = adata[:, ~adata.var_names.str.match(settings.gene_exclude_pattern, case=False, na=False)].copy()
    return adata


# Read input SpatialData
sdata = sd.read_zarr(par["input"])
sopa.make_image_patches(sdata, patch_width=par["patch_width"], patch_overlap=par["patch_overlap"])

transcript_patch_args = {
    "sdata": sdata,
    "write_cells_centroids": True,
    "patch_width": par["transcript_patch_width"],
}
transcript_patch_args["prior_shapes_key"] = par["shapes_key"]

sopa.make_transcript_patches(**transcript_patch_args)

config = {
    "dict_scale": {"x": 1, "y": 1, "z": 1},
    "mean_cell_diameter": par["mean_cell_diameter"],
    "max_cell_radius": par["max_cell_radius"],
    "norm_vector": par["norm_vector"],
    "alpha": par["alpha"], 
    "allow_disconnected_polygon": par["allow_disconnected_polygon"],
    "min_rna_per_cell": par["min_rna_per_cell"],
    "gene_column": par["gene_column"],
}


sopa.aggregation.transcripts._count_transcripts_aligned = fixed_count_transcripts_aligned
sopa.segmentation.comseg(sdata, config)

# Create output SpatialData 
sd_output = sd.SpatialData()

cell_id_col = sdata["transcripts"][f"cell_id"]
sdata.tables["table"]=ad.AnnData(obs=pd.DataFrame({"cell_id":cell_id_col}), var=sdata.tables["table"].var[[]])
sdata_new = sd.SpatialData(
    points=sdata.points,  
    tables=sdata.tables   
) 

output_path = par['output']
sdata_new.write(output_path, overwrite=True)


