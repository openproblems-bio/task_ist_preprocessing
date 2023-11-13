#!/usr/bin/env python

import txsim as tx
import scanpy as sc
import pandas as pd
import numpy as np
import os.path
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain assignments_.csv')
    parser.add_argument('-s', '--singlecell', required=True, type=str,
	    help='Path to the single cell anndata')
    parser.add_argument('-as', '--assignment', required=True, type=str, 
        help='Method list after assignments_')
    parser.add_argument('-n', '--normalize', default='total', type=str,
        help='Method to normalize raw count matrices by') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of method-specific parameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters') 
    # parser.add_argument('-t', '--threshold', default=None,
    #     help='Threshold for percent of spots with prior cell type to assign new cell type') 
    # parser.add_argument('-c', '--ctmethod', default='ssam', type=str,
    #     help='Cell type assignment method (ssam, majority, pciSeq)')
    # parser.add_argument('-ct', '--ctcertthresh', default='0.7', type=str,
    #     help='Cell type certainty threshold')
    # parser.add_argument('-g', '--pergenecorr', type=str, default='True',
    #     help='Run per gene correction')
    # parser.add_argument('-l', '--genecorrlayer', default='lognorm', type=str,
    #     help='Layer to do per gene correction on')

    
    args = parser.parse_args()

    assignment_method = args.assignment
    data = args.data
    normalize_by = args.normalize
    id_code = args.id_code
    file_sc = args.singlecell

    groupparams = eval(args.groupparams)
    if groupparams is None: groupparams = {}

    hyperparams = eval(args.hyperparams)
    if hyperparams is None: hyperparams = {}
    alpha = hyperparams.get('alpha') is not None
    use_max_area = hyperparams.get('use_max_area') is None or bool(hyperparams['use_max_area'])
    find_area = hyperparams.get('find_area') is not None and bool(hyperparams['find_area'])
    if hyperparams.get('alpha') is None: hyperparams['alpha'] = 0

    qc_params = {}
    if groupparams.get('min_counts') is not None: qc_params['min_counts'] = groupparams.get('min_counts')
    if groupparams.get('min_cell_percentage') is not None: qc_params['min_cell_percentage'] = groupparams.get('min_cell_percentage')
    qc_params['min_area'] = groupparams.get('min_area') # if None, will be handled later
    qc_params['max_area'] = groupparams.get('max_area')
        
    # Read in the single-cell data
    adata_sc = sc.read(file_sc)
    
    adata = tx.preprocessing.generate_adata(
    input_spots=pd.read_csv(f'{data}/assignments_{assignment_method}.csv'),
        adata_sc=adata_sc)
    
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
        if 'area' in adata.obs and use_max_area:
            adata.obs['area'] = np.maximum(adata.obs['alpha_area'], adata.obs['area'])
    
    if 'area' not in adata.obs and 'alpha_area' in adata.obs:
        adata.obs['area'] = adata.obs['alpha_area']

    #Normalize by area or by total counts
    if(normalize_by == 'area'):
        tx.preprocessing.normalize_by_area(adata)
    else:
        tx.preprocessing.normalize_total(adata)

    #Quality control step
    tx.preprocessing.filter_cells(adata, obs_key="passed_QC", **qc_params)
    print("Ran quality control")

    #Save AnnData object
    adata.write_h5ad(f"{data}/normcounts_{assignment_method}_{normalize_by}-{id_code}.h5ad")
