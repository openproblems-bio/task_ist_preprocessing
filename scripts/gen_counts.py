#!/usr/bin/env python

import anndata as ad
from anndata import AnnData
import pandas as pd
import scanpy as sc
import txsim as tx
import numpy as np
import os.path
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain assignments_.csv')
    parser.add_argument('-ar', '--area', default=None, type=str,
        help='Method list after areas_')
    parser.add_argument('-as', '--assignment', required=True, type=str, 
        help='Method list after assignments_')
    parser.add_argument('-n', '--normalize', default='total', type=str,
        help='Method to normalize raw count matrices by') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of parameters') 
         
    
    args = parser.parse_args()

    assignment_method = args.assignment
    data = args.data
    area_method = args.area
    normalize_by = args.normalize
    id_code = args.id_code
    hyperparams = eval(args.hyperparams)
    alpha = hyperparams is not None and hyperparams.get('alpha') is not None
    max_area = hyperparams is not None and hyperparams.get('max') is not None and hyperparams['max']
    #Read assignments
    spots = pd.read_csv(f'{data}/assignments_{assignment_method}.csv')
    spots = spots[spots['cell'] != 0]

    #Generate blank, labelled count matrix
    X = np.zeros([ len(pd.unique(spots['cell'])), len(pd.unique(spots['gene'])) ])
    adata = ad.AnnData(X, dtype = X.dtype)
    adata.obs['cell_id'] = pd.unique(spots['cell'])
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = pd.unique(spots['gene'])

    #Populate matrix using assignments
    for gene in adata.var_names:
        cts = spots[spots['gene'] == gene ]['cell'].value_counts()
        adata[:, gene] = cts.reindex(adata.obs['cell_id'], fill_value = 0)

    #TODO have different index for denovo types
    #Look for cell_types
    if(os.path.exists(f'{data}/celltypes_{assignment_method}.csv')):
        temp = pd.read_csv(f'{data}/celltypes_{assignment_method}.csv', header=None, index_col = 0)
        adata.obs['cell_type'] = pd.Categorical(temp[1][adata.obs['cell_id']])

    #Find area for normalization
    found_area = False
    if area_method is not None and normalize_by == 'area':
        #Load area data from methods provided
        temp = pd.read_csv(f'{data}/areas_{area_method}.csv', header=None, index_col = 0)
        adata.obs['area'] = temp[1][adata.obs['cell_id']]
    elif normalize_by == 'area':
        #If no provided area, search through possible areas
        methods = assignment_method
        method_list = assignment_method.split('_')
        #Work backwards through method list until areas file is found
        for i in range(0, len(method_list)):
            methods = '_'.join(method_list)
            if(os.path.exists(f'{data}/areas_{methods}.csv')):
                temp = pd.read_csv(f'{data}/areas_{methods}.csv', header=None, index_col = 0)
                adata.obs['area'] = temp[1][adata.obs['cell_id']]
                found_area = True
                break
            method_list.pop()
        if not found_area:
            #If none found, use alpha area
            alpha = True        
    
    # Calculate area based on alpha shape from molecules for each shape
    # If there are <3 molecules for a cell, use the mean area per molecule
    # times the number of molecules in the cell
    if alpha and normalize_by == 'area':
        import alphashape
        from descartes import PolygonPatch
        area_vec = np.zeros([adata.n_obs])
        for i in range(adata.n_obs):
            dots = pd.concat(
                [spots[spots['cell'] == adata.obs['cell_id'][i]].x,
                spots[spots['cell'] == adata.obs['cell_id'][i]].y],
                axis=1
            )
            pts = list(dots.itertuples(index=False, name=None))
            if(len(pts) > 2):
                alpha_shape = alphashape.alphashape(pts,hyperparams['alpha'])
                area_vec[i] = alpha_shape.area
            else:
                area_vec[i] = np.nan
        #Normalize each cell by the number of molecules assigned to it
        mean_area = np.nanmean(area_vec / np.sum(adata.X, axis=1) )
        #Use this mean area to fill in NaN values
        area_vec[np.isnan(area_vec)] = mean_area * np.sum(adata.X, axis=1)[np.isnan(area_vec)]
        adata.obs['alpha_area'] = area_vec
        if area_method is not None and max_area:
            adata.obs['area'] = np.maximum(adata.obs['alpha_area'], adata.obs['area'])
        elif area_method is None:
            adata.obs['area'] = adata.obs['alpha_area']
    
    #Normalize by area or by total counts
    if(normalize_by == 'area'):
        tx.preprocessing.normalize_by_area(adata)
    else:
        tx.preprocessing.normalize_total(adata)

    #Save AnnData object
    adata.write_h5ad(f"{data}/counts_{assignment_method}_{normalize_by}-{id_code}.h5ad")
