import spatialdata as sd
from tifffile import imwrite  
import sys
import numpy as np
import pandas as pd
import xarray as xr
import dask
import txsim as tx
import anndata as ad
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process input files and generate output files for FastReseg input"
    )
    parser.add_argument('input_path_ist', help='Path to the input file')
    parser.add_argument('input_segmentation_path', help='Path to the input segmentation file')
    parser.add_argument('input_sc_reference_path', help='Path to the input single-cell reference file')
    parser.add_argument('output_path_counts', help='Path for the output TSV file')
    parser.add_argument('output_path_transcripts', help='Path for the output TIF file')
    parser.add_argument('output_path_cell_type', help='Output cell type specification')
    
    return parser.parse_args()

### parsing arguments
args = parse_arguments()
print("args:")
print(args)
input_path = args.input_path_ist
input_segmentation_path = args.input_segmentation_path
input_sc_reference_path = args.input_sc_reference_path
print("path")
print(input_sc_reference_path)
output_path_counts = args.output_path_counts
output_path_transcripts = args.output_path_transcripts
output_path_cell_type = args.output_path_cell_type

## potential other parameters (TODO - make configurable)
um_per_pixel = 0.5
sc_celltype_key = 'cell_type'

### reading the data in
sdata = sd.read_zarr(input_path)

### reading in basic segmentation
sdata_segm = sd.read_zarr(input_segmentation_path)
segmentation_coord_systems = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).keys()

# In case of a translation transformation of the segmentation (e.g. crop of the data), we need to adjust the transcript coordinates
trans = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True)['global'].inverse()

transcripts = sd.transform(sdata['transcripts'], to_coordinate_system='global')
transcripts = sd.transform(transcripts, trans, 'global')

print('Assigning transcripts to cell ids', flush=True)
y_coords = transcripts.y.compute().to_numpy(dtype=np.int64)
x_coords = transcripts.x.compute().to_numpy(dtype=np.int64)
if isinstance(sdata_segm["segmentation"], xr.DataTree):
    label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy() 
else:
    label_image = sdata_segm["segmentation"].to_numpy()
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        label_image[y_coords, x_coords], chunks=tuple(sdata['transcripts'].map_partitions(len).compute())
    ), 
    index=sdata['transcripts'].index
)
sdata['transcripts']["cell_id"] = cell_id_dask_series

### extracting transcript ids
print('Transforming transcripts coordinates', flush=True)
transcripts = sd.transform(sdata['transcripts'], to_coordinate_system='global')

transcripts_df = transcripts.compute()
transcripts_df.rename(columns = {'feature_name': 'target', 
'transcript_id': 'UMI_transID', 'cell_id': 'UMI_cellID'}, inplace = True)

transcripts_df = transcripts_df.loc[:, ['target', 'x', 'y', 'z', 'UMI_transID', 'UMI_cellID']]
transcripts_df.to_csv(output_path_transcripts)


#### aggregating counts per transcript, based on 
df = sdata['transcripts'].compute()
df.feature_name = df.feature_name.astype(str)

adata_sp = tx.preprocessing.generate_adata(df, cell_id_col='cell_id', gene_col='feature_name') #TODO: x and y refers to a specific coordinate system. Decide which space we want to use here. (probably should be handled in the previous assignment step)
adata_sp.layers['counts'] = adata_sp.layers['raw_counts']
del adata_sp.layers['raw_counts']
adata_sp.var["gene_name"] = adata_sp.var_names
print(adata_sp.var_names[1:10])

# currently the function also saves the transcripts in the adata object, but this is not necessary here
del adata_sp.uns['spots']
del adata_sp.uns['pct_noise']


count_df = pd.DataFrame(adata_sp.X.toarray(), 
                       index=adata_sp.obs_names, 
                       columns=adata_sp.var_names)
count_df.to_csv(output_path_counts)

#### run cell annotation with ssam
adata_sc = ad.read_h5ad(input_sc_reference_path)
adata_sc.X = adata_sc.layers["normalized"]
print(adata_sc.var_names[1:10])

shared_genes = [g for g in adata_sc.var_names if g in adata_sp.var_names]
adata_sp = adata_sp[:,shared_genes] 

print('Annotating cell types', flush=True)
adata_sp = tx.preprocessing.run_ssam(
    adata_sp, transcripts.compute(), adata_sc, um_p_px=um_per_pixel, 
    cell_id_col='cell_id', gene_col='feature_name', sc_ct_key=sc_celltype_key
)
cell_type_df = adata_sp.obs["ct_ssam"].astype(str)
cell_type_df.to_csv(output_path_cell_type, header=True)