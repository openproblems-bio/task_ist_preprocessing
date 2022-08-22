#!/usr/bin/env python

import txsim as tx
import pandas as pd
import numpy as np
import os.path
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain assignments_.csv')
    parser.add_argument('-as', '--assignment', required=True, type=str, 
        help='Method list after assignments_')
    parser.add_argument('-n', '--normalize', default='total', type=str,
        help='Method to normalize raw count matrices by') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of parameters') 
    parser.add_argument('-t', '--threshold', default=None,
        help='Threshold for percent of spots with prior cell type to assign new cell type') 
    
    args = parser.parse_args()

    assignment_method = args.assignment
    data = args.data
    normalize_by = args.normalize
    id_code = args.id_code
    hyperparams = eval(args.hyperparams)
    if hyperparams is None: hyperparams = {}
    alpha = hyperparams.get('alpha') is not None
    max_area = hyperparams.get('max') is None or hyperparams['max']
    find_area = hyperparams.get('find_area') is not None and hyperparams['find_area']
    prior_pct = 0.7 if args.threshold is None else args.threshold 
    if hyperparams.get('alpha') is None: hyperparams['alpha'] = 0

    #TODO have different index for denovo types
    #Look for cell_types
    if os.path.exists(f'{data}/celltypes_{assignment_method}.csv'):
        adata = tx.preprocessing.generate_adata(
            molecules=f'{data}/assignments_{assignment_method}.csv',
            cell_types=f'{data}/celltypes_{assignment_method}.csv',
            prior_pct = prior_pct)
    else:
        adata = tx.preprocessing.generate_adata(
            molecules=f'{data}/assignments_{assignment_method}.csv', 
            prior_pct=prior_pct)

    #Find area for normalization
    if normalize_by == 'area' or find_area:
        methods = assignment_method
        method_list = assignment_method.split('_')
        #Work backwards through method list until areas file is found
        for i in range(0, len(method_list)):
            methods = '_'.join(method_list)
            if(os.path.exists(f'{data}/areas_{methods}.csv')):
                temp = pd.read_csv(f'{data}/areas_{methods}.csv', header=None, index_col = 0)
                adata.obs['area'] = temp[1][adata.obs['cell_id']].values
                break
            method_list.pop()

    #If none found, use alpha area
    if 'area' not in adata.obs:
        alpha = True        
    
    # Calculate area based on alpha shape from molecules for each shape
    if alpha and (normalize_by == 'area' or find_area):
        tx.preprocessing.calculate_alpha_area(
            adata=adata,
            alpha=hyperparams['alpha']
        )
        if 'area' in adata.obs and max_area:
            adata.obs['area'] = np.maximum(adata.obs['alpha_area'], adata.obs['area'])
    
    if 'area' not in adata.obs and 'alpha_area' in adata.obs:
        adata.obs['area'] = adata.obs['alpha_area']

    #Normalize by area or by total counts
    if(normalize_by == 'area'):
        tx.preprocessing.normalize_by_area(adata)
    else:
        tx.preprocessing.normalize_total(adata)


    #Save AnnData object
    adata.write_h5ad(f"{data}/counts_{assignment_method}_{normalize_by}-{id_code}.h5ad")
