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
    parser.add_argument('-c', '--counts', required=True, type=str, 
        help='Method list after counts_')
    parser.add_argument('-a', '--annotate', required=True, type=str,
        help='Method to annotate celltypes') 
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

    counts_method = args.counts
    data = args.data
    annotate_with = args.annotate
    id_code = args.id_code
    file_sc = args.singlecell

    groupparams = eval(args.groupparams)
    hyperparams = eval(args.hyperparams)
    if groupparams is None: groupparams = {}
    if hyperparams is None: hyperparams = {}

    per_gene_correction = groupparams.get('per_gene_correction') is not None and bool(groupparams['per_gene_correction'])
    gene_corr_layer = groupparams.get('gen_corr_layer') if groupparams.get('gen_corr_layer') is not None else 'lognorm'
    prior_celltypes = None

    # Read in the single-cell data
    adata_sc = sc.read(file_sc)
    
    adata = sc.read(f'{data}/counts_{counts_method}.h5ad')

    if annotate_with == 'pciSeqCT':
        methods = counts_method
        method_list = counts_method.split('_')
        #Work backwards through method list until areas file is found
        for i in range(0, len(method_list)):
            methods = '_'.join(method_list)
            if(os.path.exists(f'{data}/celltypes_{methods}.csv')):
                prior_celltypes = pd.read_csv(f'{data}/celltypes_{methods}.csv', index_col = 0)
                break
            method_list.pop()

    adata = tx.preprocessing.annotate_celltypes(
        adata,
        adata_sc = adata_sc,
        ct_method=annotate_with,
        hyperparams=hyperparams,
        prior_celltypes=prior_celltypes,
    )

    # Do per-gene correction if active
    
    if per_gene_correction:
        tx.preprocessing.gene_efficiency_correction(adata, adata_sc, gene_corr_layer)
        if gene_corr_layer!='lognorm':
            adata.layers['lognorm'] = adata.layers['norm']
            sc.pp.log1p(adata, layer='lognorm')


    #Save AnnData object
    adata.write_h5ad(f"{data}/counts_{counts_method}_{annotate_with}-{id_code}.h5ad")
