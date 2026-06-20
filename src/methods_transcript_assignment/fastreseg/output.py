import spatialdata as sd
from tifffile import imwrite
import sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.sparse import coo_matrix
from pathlib import Path
import os
import dask
import anndata as ad

def convert_to_lower_dtype(arr):
    max_val = arr.max()
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64

    return arr.astype(new_dtype)
    
path_to_cell_ids_out = sys.argv[1]
path_to_gene_names_out = sys.argv[2] 
path_to_transcripts_out = sys.argv[3]

output_zarr_path = sys.argv[4]

cell_df = pd.read_csv(path_to_cell_ids_out, index_col=0)
var_names = pd.read_csv(path_to_gene_names_out, index_col=0)
df = pd.read_csv(path_to_transcripts_out, index_col=0)
print(df.head())

##converting to dask transcript output for spatialdata format
df = df.loc[:,["x", "y", "z", "target", "updated_cellID", "updated_celltype", "UMI_transID"]]
df.rename(columns={"updated_cellID": "cell_id", "target": "feature_name", "UMI_transID": "transcript_id"}, inplace=True)
transcripts_dask = dask.dataframe.from_pandas(df, npartitions = 1)


sdata_transcripts_only = sd.SpatialData(
    points={
        "transcripts": sd.models.PointsModel.parse(transcripts_dask)
    },
    tables={
        "table": ad.AnnData(
          obs=cell_df.loc[:,['updated_cellID', 'updated_celltype']],
          var=pd.DataFrame(index=var_names['x'])
        )
    }
)
sdata_transcripts_only.write(output_zarr_path)