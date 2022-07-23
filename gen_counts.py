#!/usr/bin/env python

import anndata as ad
from anndata import AnnData
import pandas as pd
import scanpy as sc
import txsim as tx
import numpy as np
import os.path
import argparse

#INPUT: assignments.csv, [optional] area.csv
#OUTPUT: counts.h5ad
#From config:
#segmentation_method = 'imagej'
#assignment_method = 'pciSeq'
#area_method = 'alpha' 
#normalize_by = 'area'
#TODO add method for taking max of several methods

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-s', '--segment', required=True, type=str,
        help='Segmentation method used for image')
    parser.add_argument('-a', '--assignment', required=True, type=str, 
        help='Assignment method used for molecules')
    parser.add_argument('-n', '--normalize', default='total', type=str,
        help='Method to normalize raw count matrices by') 
    parser.add_argument('-am', '--areamethod', default=None, type=str,
        help='Method to calculate area for each cell, leave as None to use previous areas') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
         
    
    args = parser.parse_args()

    assignment_method = args.assignment
    data = args.data
    segmentation_method = args.segment
    normalize_by = args.normalize
    area_method = args.areamethod
    id_code = args.id_code

    #Read assignments
    spots = pd.read_csv(f'{data}/assignments_{segmentation_method}_{assignment_method}.csv')
    spots = spots[spots['cell'] != 0]

    #Generate blank, labelled count matrix
    X = np.empty([ len(pd.unique(spots['cell'])), len(pd.unique(spots['gene'])) ])
    adata = ad.AnnData(X, dtype = X.dtype)
    adata.obs_names = pd.unique(spots['cell'])
    adata.var_names = pd.unique(spots['gene'])

    #Populate matrix using assignments
    for index, row in spots.iterrows():
        n = row['cell']
        g = row['gene']
        adata[adata.obs_names==n, g] = adata[adata.obs_names==n, g].to_df()[g][n] + 1


    #Load area data from same as assignment or segmentation method
    if(area_method is None):
        if(os.path.exists(f'{data}/areas_{assignment_method}.csv')):
            temp = pd.read_csv(f'{data}/areas_{assignment_method}.csv', header=None)
            adata.obs['area'] = temp[1][adata.obs_names]
        elif(os.path.exists(f'{data}/areas_{segmentation_method}.csv')):
            temp = pd.read_csv(f'{data}/areas_{segmentation_method}.csv', header=None)
            adata.obs['area'] = temp[1][adata.obs_names]
        else:
            #If no area data detected, use alpha area from points
            area_method = 'alpha'

    # Calculate area based on alpha shape from molecules for each shape
    # If there are <3 molecules for a cell, use the mean area per molecule
    # times the number of molecules in the cell
    if(area_method == 'alpha'):
        import alphashape
        from descartes import PolygonPatch
        area_vec = np.empty([adata.n_obs])
        for i in range(adata.n_obs):
            dots = pd.concat(
                [spots[spots['cell'] == adata.obs_names[i]].x,
                spots[spots['cell'] == adata.obs_names[i]].y],
                axis=1
            )
            pts = list(dots.itertuples(index=False, name=None))
            if(len(pts) > 2):
                alpha_shape = alphashape.alphashape(pts,0.)
                area_vec[i] = alpha_shape.area
            else:
                area_vec[i] = np.nan
        #Normalize each cell by the number of molecules assigned to it
        mean_area = np.nanmean(area_vec / np.sum(adata.X, axis=1) )
        #Use this mean area to fill in NaN values
        area_vec[np.isnan(area_vec)] = mean_area * np.sum(adata.X, axis=1)[np.isnan(area_vec)]
        adata.obs['area'] = area_vec
    elif (area_method is not None):
        temp = pd.read_csv(f'{data}/areas_{area_method}.csv', header=None)
        adata.obs['area'] = temp[1][adata.obs_names]
    
    #Normalize by area or by total counts
    if(normalize_by == 'area'):
        tx.preprocessing.normalize_by_area(adata)
    else:
        tx.preprocessing.normalize_total(adata)

    #Save AnnData object
    adata.write_h5ad(f"{data}/counts_{segmentation_method}_{assignment_method}_{normalize_by}-{id_code}.h5ad")
    print(f'Saved {data}/counts_{segmentation_method}_{assignment_method}_{normalize_by}-{id_code}.h5ad')
