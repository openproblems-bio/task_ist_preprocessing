import spatialdata as sd
from tifffile import imwrite
import sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.sparse import coo_matrix
from pathlib import Path
import os

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
    
input_path = sys.argv[1]
input_cellspot = sys.argv[2]
output_path = sys.argv[3]

print(input_path)

sdata = sd.read_zarr(input_path)
sd_output = sd.SpatialData()

transformation = sdata['morphology_mip']['scale0'].image.transform.copy()

scs_output = pd.read_csv(input_cellspot, sep = "\t", header = None)
scs_output.rename(columns = {0: 'spot', 1: 'cell'}, inplace = True)
scs_output[['x', 'y']] = scs_output.spot.str.split(':', n=1, expand=True)

scs_output.x = scs_output.x.astype(int)
scs_output.y = scs_output.y.astype(int)

scs_output.cell = scs_output.cell.astype(int)

##converting back to image-like format
sparse_matrix = coo_matrix(
        (scs_output['cell'], (scs_output['x'], scs_output['y'])), 
        shape=(scs_output.x.max() + 1, scs_output.y.max() + 1)
    )
labels = sparse_matrix.toarray()

## all non-segmented as -1
labels = np.where(labels == 0, -1, labels)
labels = convert_to_lower_dtype(labels)

labels_array =xr.DataArray(labels, name=f'segmentation', dims=('y', 'x'))
parsed_labels = sd.models.Labels2DModel.parse(labels_array, transformations=transformation)
sd_output.labels['segmentation'] = parsed_labels

print("Writing output", flush=True)
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
if os.path.exists(output_path):
    shutil.rmtree(output_path)
sd_output.write(output_path)