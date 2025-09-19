import dask
import xarray as xr
import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import numpy as np

## VIASH START
par = {
    "input_ist": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
    "input_segmentation": "resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr",
    "transcripts_key": "transcripts",
    "coordinate_system": "global",
    "output": "temp/comseg/transcripts.zarr",

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



# Read input files
print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])


# Convert the prior segmentation to polygons
if isinstance(sdata_segm["segmentation"], xr.DataTree):
    shapes_gdf = sopa.shapes.vectorize(sdata_segm["segmentation"]["scale0"].image)
else:
    shapes_gdf = sopa.shapes.vectorize(sdata_segm["segmentation"])

sdata["segmentation_boundaries"] = sd.models.ShapesModel.parse(
    shapes_gdf, transformations=sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).copy()
)

# Make patches
sopa.make_image_patches(sdata, patch_width=par["patch_width"], patch_overlap=par["patch_overlap"])

transcript_patch_args = {
    "sdata": sdata,
    "write_cells_centroids": True,
    "patch_width": par["transcript_patch_width"],
    "prior_shapes_key": "segmentation_boundaries",
}

sopa.make_transcript_patches(**transcript_patch_args)

# Run ComSeg
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
# sopa.settings.parallelization_backend = 'dask'
sopa.segmentation.comseg(sdata, config)

# Assign transcripts to cell ids
sopa.spatial.assign_transcript_to_cell(
    sdata,
    points_key="transcripts",
    shapes_key="comseg_boundaries",
    key_added="cell_id",
    unassigned_value=0
)

# Create output SpatialData 

# Create objects for cells table
print('Creating objects for cells table', flush=True)
unique_cells = np.unique(sdata["transcripts"]["cell_id"])
zero_idx = np.where(unique_cells == 0)
if len(zero_idx[0]): 
    unique_cells=np.delete(unique_cells, zero_idx[0][0])
cell_id_col = pd.Series(unique_cells, name='cell_id', index=unique_cells)

# Create transcripts only sdata
print('Subsetting to transcripts cell id data', flush=True)
sdata_transcripts_only = sd.SpatialData(
    points={
        "transcripts": sdata['transcripts']
    },
    tables={
        "table": ad.AnnData(
          obs=pd.DataFrame(cell_id_col),
          var=sdata.tables["table"].var[[]]
        )
    }
)


output_path = par['output']
sdata_transcripts_only.write(output_path, overwrite=True)


