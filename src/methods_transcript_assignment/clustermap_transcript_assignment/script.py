import numpy as np
import dask
import spatialdata as sd
import anndata as ad
import pandas as pd
import os
import shutil

import time
from datetime import timedelta
import math
from ClusterMap.clustermap import ClusterMap
from ClusterMap.utils import get_img, split


## VIASH START
par = {
  'input_ist': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
  'input_segmentation': 'resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr',
  'transcripts_key': 'transcripts',
  'coordinate_system': 'global',
  'output': 'basic_assigned_transcripts.zarr',
  'window_size': 700,
  'use_dapi': True,
  'xy_radius': 40,
  'z_radius': 0,
  'fast_preprocess': False ,
  'gauss_blur': True,
  'sigma': 1,
  'pct_filter': 0.0,
  'LOF': False,
  'contamination': 0,
  'min_spot_per_cell': 5,
  'cell_num_threshold': 0.1 ,
  'add_dapi': True,
  'use_genedis': True,
  'dapi_grid_interval': 5
}
meta = {
  'name': 'clustermap'
}
## VIASH END

# ----------------- FUNCTIONS ----------------- #

def format_spots_for_clustermap(spots):
    """
    """
    # Rename columns
    spots.rename(columns={"feature_name": 'gene_name',
                          "x": 'spot_location_1',
                          "y": 'spot_location_2',
                          "z": 'spot_location_3',
                          }, inplace=True)
    
    # Make sure that index goes from 0 to n-1
    spots.index = [i for i in range(len(spots))]
    
    # Convert gene_name to gene identity
    genes = pd.DataFrame(spots['gene_name'].unique())
    a1 = list(genes[0])
    gene = list(map(lambda x: a1.index(x)+1, spots['gene_name']))
    spots = spots.astype({'spot_location_1': int, 'spot_location_2': int})
    if 'spot_location_3' in spots.columns:
        spots = spots.astype({'spot_location_3': int})
    spots['gene'] = gene
    spots['gene'] = spots['gene'].astype('int')
    
    return spots

def run_clustermap_over_chunks(
        dapi, spots, window_size=1500, xy_radius=40, z_radius=0, gauss_blur=True, sigma=1, fast_preprocess=False, 
        dapi_grid_interval=5, pct_filter=0.0, LOF=False, contamination=0.0, min_spot_per_cell=5, 
        cell_num_threshold=0.1, add_dapi=True, use_genedis=True, 
    ):
    """
    """

    # Infer parameters
    num_gene=np.max(spots['gene'])
    num_dims = len(dapi.shape)
    gene_list=np.arange(1,num_gene+1)
    
    # Initialize model over full image
    model = ClusterMap(
        spots=spots, dapi=dapi, gene_list=gene_list, num_dims=num_dims, xy_radius=xy_radius, z_radius=z_radius, 
        fast_preprocess=fast_preprocess, gauss_blur=gauss_blur, sigma=sigma
    )
    
    # Set all spots to background
    model.spots['clustermap'] = -1
    
    # Trim to tiles
    label_img = get_img(dapi, spots, window_size=window_size, margin=math.ceil(window_size*0.1))
    tiles_df = split(dapi, label_img, spots, window_size=window_size, margin=math.ceil(window_size*0.1))
    
    # Process each tile
    print('####### Start tile 1')
    t0 = time.time()
    for tile_num, (_, _, dapi_tile, spots_tile, _) in tiles_df.iterrows():
        
        if spots_tile.shape[0] < 20:
            continue
    
        # Instantiate model for tile
        t0_model = time.time()
        model_tile = ClusterMap(
            spots=spots_tile, dapi=dapi_tile, gene_list=gene_list, num_dims=num_dims, xy_radius=xy_radius, 
            z_radius=z_radius, fast_preprocess=fast_preprocess, gauss_blur=gauss_blur,sigma=sigma
        )
        time_model = timedelta(seconds=(time.time() - t0_model))
        
    
        # Preprocessing
        t0_prepro = time.time()
        model_tile.preprocess(
            dapi_grid_interval=dapi_grid_interval, pct_filter=pct_filter, LOF=LOF, contamination=contamination
        )
        time_prepro = timedelta(seconds=(time.time() - t0_prepro))
    
        # Segmentation
        model_tile.min_spot_per_cell=min_spot_per_cell
        t0_segment = time.time()
        model_tile.segmentation(
            cell_num_threshold=cell_num_threshold, dapi_grid_interval=dapi_grid_interval, add_dapi=add_dapi, 
            use_genedis=use_genedis
        )
        time_segment = timedelta(seconds=(time.time() - t0_segment))
        
        # Stitch new tile to the previous ones
        t0_stitch = time.time()
        model.stitch(model_tile, tiles_df, tile_num)
        time_stitch = timedelta(seconds=(time.time() - t0_stitch))
        
        
        time_formatted = timedelta(seconds=(time.time() - t0))
        expected_time = timedelta(seconds=(time.time() - t0) / (tile_num+1) * len(tiles_df)) if tile_num > 0 else "N/A"
        print(
            f'####### tile: {tile_num+1}/{len(tiles_df)} | time: {time_formatted} | expected total time: {expected_time}'+
            f'\n\t time init model: {time_model} \n\t time preprocess: {time_prepro} ' + 
            f'\n\t time segmentation: {time_segment} \n\t time stitch: {time_stitch}',
            flush=True
        )
        #TODO: include an assertion for an expected_time < 24 h? And a recommendation to decrease dapi_grid_interval?/data_size
        
    
    return model.spots

def run_clustermap(
        dapi, spots, xy_radius=40, z_radius=0, gauss_blur=True, sigma=1, fast_preprocess=False, dapi_grid_interval=5, 
        pct_filter=0.0, LOF=False, contamination=0.0, min_spot_per_cell=5, cell_num_threshold=0.1, add_dapi=True,  
        use_genedis=True,
    ):
    """"""
    
    # Infer parameters
    num_gene = np.max(spots['gene'])
    num_dims = len(dapi.shape)
    gene_list = np.arange(1, num_gene+1)
    
    # Create Model
    model = ClusterMap(
        spots=spots, dapi=dapi, gene_list=gene_list, num_dims=num_dims, xy_radius=xy_radius, z_radius=z_radius, 
        fast_preprocess=fast_preprocess, gauss_blur=gauss_blur, sigma=sigma
    )

    # Preprocessing
    model.preprocess(dapi_grid_interval=dapi_grid_interval, pct_filter=pct_filter, LOF=LOF, contamination=contamination)
    
    # Segmentation
    model.min_spot_per_cell = min_spot_per_cell
    model.segmentation(
        cell_num_threshold=cell_num_threshold, dapi_grid_interval=dapi_grid_interval, add_dapi=add_dapi, 
        use_genedis=use_genedis
    )

    return model.spots



# ----------------- SCRIPT ----------------- #



assert par['input_segmentation'], 'Segmentation input is required for this assignment method.'

print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])

# Check if coordinate system is available in input data
transcripts_coord_systems = sd.transformations.get_transformation(sdata[par["transcripts_key"]], get_all=True).keys()
assert par['coordinate_system'] in transcripts_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."
segmentation_coord_systems = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).keys()
assert par['coordinate_system'] in segmentation_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."

print('Transforming transcripts coordinates', flush=True)
transcripts = sd.transform(sdata[par['transcripts_key']], to_coordinate_system=par['coordinate_system'])

# In case of a translation transformation of the segmentation (e.g. crop of the data), we need to adjust the transcript coordinates
trans = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True)[par['coordinate_system']].inverse()
transcripts = sd.transform(transcripts, trans, par['coordinate_system'])

print('Assigning transcripts to cell ids', flush=True)

# Turn segmentation SpatialData into numpy array
label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy() #TODO: mabye this line needs generalization (DataTree vs DataArray)
dapi_image = np.squeeze(sdata['morphology_mip']['scale0']['image'].compute())

# Extract coordinates and feature names (= gene names) from SpatialData
# and convert into the ClusterMap format
spots_cmap = format_spots_for_clustermap(transcripts.compute()[['feature_name', 'x', 'y', 'z']])

# Run clustermap
#TODO: there is probably a more effecient way to assign the transcripts without using so many dataframes
# however this is copied over from the old pipeline

#TODO CLUSTERMAP NEEDS RAW DAPI NOT LABELLED IMAGE

if par["window_size"] is None:
    spots_cmap = run_clustermap(dapi_image, spots_cmap,
        xy_radius=par["xy_radius"], 
        z_radius=par["z_radius"], 
        gauss_blur=par["gauss_blur"], 
        sigma=par["sigma"], 
        fast_preprocess=par["fast_preprocess"], 
        dapi_grid_interval=par["dapi_grid_interval"], 
        pct_filter=par["pct_filter"], 
        LOF=par["LOF"], 
        contamination=par["contamination"], 
        min_spot_per_cell=par["min_spot_per_cell"], 
        cell_num_threshold=par["cell_num_threshold"], 
        add_dapi=par["xy_radius"],  
        use_genedis=par["xy_radius"])
else:
    spots_cmap = run_clustermap_over_chunks(dapi_image, spots_cmap,
        window_size=par["window_size"],
        xy_radius=par["xy_radius"], 
        z_radius=par["z_radius"], 
        gauss_blur=par["gauss_blur"], 
        sigma=par["sigma"], 
        fast_preprocess=par["fast_preprocess"], 
        dapi_grid_interval=par["dapi_grid_interval"], 
        pct_filter=par["pct_filter"], 
        LOF=par["LOF"], 
        contamination=par["contamination"], 
        min_spot_per_cell=par["min_spot_per_cell"], 
        cell_num_threshold=par["cell_num_threshold"], 
        add_dapi=par["xy_radius"],  
        use_genedis=par["xy_radius"])

# Correct cell_ids since ClusterMap sets -1 instead of 0 as noise
spots_cmap['clustermap'] += 1 

#assign ClusterMap output (stored in spots_cmap) to the transcripts (i.e. assign transcripts to cells)
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        spots_cmap["clustermap"].values, chunks=tuple(sdata[par['transcripts_key']].map_partitions(len).compute())
    ), 
    index=sdata[par['transcripts_key']].index
)
sdata[par['transcripts_key']]["cell_id"] = cell_id_dask_series 

#create new .obs for cells based on the segmentation output (corresponding with the transcripts 'cell_id')
unique_cells = np.unique(spots_cmap["clustermap"].values)

# check if a '0' (noise/background) cell is in cell_id and remove
zero_idx = np.where(unique_cells == 0)
if len(zero_idx[0]): unique_cells=np.delete(unique_cells, zero_idx[0][0])

#transform into pandas series and check
cell_id_col = pd.Series(unique_cells, name='cell_id', index=unique_cells)
assert 0 not in cell_id_col, "Found '0' in cell_id column of assingment output cell matrix"


# TODO: Also take care of the following cases:
# - segmentation 3D, transcripts 3D
# - segmentation 3D, transcripts 2D
# - segmentation 2D, transcripts 3D

print('Subsetting to transcripts cell id data', flush=True)
sdata_transcripts_only = sd.SpatialData(
  points={
    "transcripts": sdata[par['transcripts_key']]
  },
  tables={
    "table": ad.AnnData(
      obs=pd.DataFrame(cell_id_col),
      var=sdata.tables["table"].var[[]]
    )
  }
)

print('Write transcripts with cell ids', flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_transcripts_only.write(par['output'])
