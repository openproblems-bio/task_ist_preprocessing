#!/usr/bin/env python

import txsim as tx
import scanpy as sc
import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-sc', '--singlecell', required=True, type=str, 
        help='Non-normalized adata')
    #parser.add_argument('-m', '--molecules', required=True, type=str, 
    #    help='Spatial spots csv')
    parser.add_argument('-m', '--molecules', required=True, nargs='+', type=str,
                    help='Spatial spots csv (expects a list of strings)')
    parser.add_argument('-o', '--output', default='total', type=str,
        help='Adata output') 
    
    args = parser.parse_args()

    file_in = args.singlecell
    file_spots_list = args.molecules
    file_out = args.output

    adata = sc.read(file_in)
    adata = tx.preprocessing.normalize_sc(adata)
    
    # NOTE: Arbitrary columns can lead to issues for anndata in R (specifically in the method script for 
    #       scrattch.mapping). In the future we'll need more columns (region matching filter) --> TODO
    adata.obs = adata.obs[["celltype"]]
    
    genes_sp = []
    for file_spots in file_spots_list:
        spots = pd.read_csv(file_spots)
        genes_sp = list(np.unique(genes_sp + spots["Gene"].unique().tolist()))
        
    adata = adata[:,[g for g in adata.var_names if g in genes_sp]]
    
    adata.write_h5ad(file_out)
