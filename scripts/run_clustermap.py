#!/usr/bin/env python

from pathlib import Path
import argparse
import yaml
import tifffile
import time
from datetime import timedelta
import math
import numpy as np
import pandas as pd
from ClusterMap.clustermap import ClusterMap
from ClusterMap.utils import get_img, split

# ----------------- FUNCTIONS ----------------- #

def format_spots_for_clustermap(spots):
    """
    """
    # Rename columns
    spots.rename(columns={"Gene": 'gene_name',
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells using ClusterMap')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    # Load default hyperparameter
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        hparams_defaults = yaml.safe_load(file)["clustermap"]
    
    # Get input arguments
    args = parser.parse_args()
    molecules = args.molecules
    data = args.data
    image = args.input
    id_code = args.id_code
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    print("###################\n", hyperparams, flush=True)
    
    #Read input data
    dapi = tifffile.imread(image)
    spots = pd.read_csv(molecules)

    # Format spots for clustermap (spots_cm) and keep a copy (spots)
    spots_cm = spots.copy()
    spots_cm = format_spots_for_clustermap(spots_cm)
    
    # Run clustermap
    if hyperparams["window_size"] is None:
        hyperparams.pop("window_size")
        spots_cm = run_clustermap(dapi, spots_cm, **hyperparams)
    else:
        spots_cm = run_clustermap_over_chunks(dapi, spots_cm, **hyperparams)
    
    # Add clustermap cell assignment results to spots
    spots["cell"] = spots_cm["clustermap"].values
    spots["cell"] += 1
    
    #Save assignments to csv
    spots.to_csv(f'{data}/assignments_clustermap-{id_code}.csv', index = False)
