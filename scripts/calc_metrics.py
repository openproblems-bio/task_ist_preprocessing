#!/usr/bin/env python
import pandas as pd
import numpy as np
import scanpy as sc
import txsim as tx
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also count matrix')
    parser.add_argument('-m', '--methods', required=True, type=str,
        help='Methods used to generate count matrix')
    parser.add_argument('-sc', '--singlecell', required=True, type=str,
        help='Single cell h5ad count matrix with celltype in anndata.obs[\'celltype\']') 
    
    args = parser.parse_args()
    
    methods = args.methods
    data = args.data
    sc_data = args.singlecell

    #Read count matrices
    scdata = sc.read(sc_data)
    stdata = sc.read(f'{data}/counts_{methods}.h5ad')

    #Gather information on spatial data
    stdata.obs['n_counts']= np.sum(stdata.layers['raw_counts'], axis = 1)
    stdata.obs['n_unique_genes']= np.sum(stdata.layers['raw_counts']>0, axis = 1)
    stdata.var['n_counts']= np.sum(stdata.layers['raw_counts'], axis=0)
    stdata.var['n_unique_cells']= np.sum(stdata.layers['raw_counts']>0, axis = 0)
    
    #Generate metrics
    metrics = {}
    metrics['coex_all'] = tx.metrics.coexpression_similarity(stdata, scdata)
    metrics['coex_thresh'] = tx.metrics.coexpression_similarity(stdata, scdata, thresh=0.5)
    ct =  tx.metrics.coexpression_similarity_celltype(stdata, scdata, thresh=0)
    metrics['coex_bytype_thresh'] = np.nanmean(ct['mean_diff'])
    idx = ~np.isnan(ct['mean_diff'])
    metrics['coex_bytype_weighted_thresh'] = np.average(ct['mean_diff'][idx], weights = ct['pct'][idx])
    metrics['pct_spots_unassigned'] = stdata.uns['pct_noise']
    metrics['n_cells'] = stdata.n_obs
    metrics['mean_cts_per_cell'] = np.mean(stdata.obs['n_counts'])
    metrics['mean_genes_per_cell'] = np.mean(stdata.obs['n_unique_genes'])
    metrics['mean_cts_per_gene'] = np.mean(stdata.var['n_counts'])
    metrics['mean_cells_per_gene'] = np.mean(stdata.var['n_unique_cells'])
    metrics['pct_cells_no_type'] = stdata.obs['celltype'].value_counts()['None'] / stdata.n_obs
    

    df = pd.DataFrame.from_dict(metrics, orient='index')

    df.to_csv(f'{data}/metrics_{methods}.csv')
    