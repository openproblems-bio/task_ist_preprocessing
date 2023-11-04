#!/usr/bin/env python

import txsim as tx
import scanpy as sc
import pandas as pd
import os.path
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-sc', '--singlecell', required=True, type=str, 
        help='Non-normalized adata')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Spatial spots csv')
    parser.add_argument('-o', '--output', default='total', type=str,
        help='Adata output') 
    
    args = parser.parse_args()

    file_in = args.singlecell
    file_spots = args.molecules
    file_out = args.output

    adata = sc.read(file_in)
    adata = tx.preprocessing.normalize_sc(adata)
    
    spots = pd.read_csv(file_spots)
    genes_sp = spots["Gene"].unique().tolist()
    adata = adata[:,[g for g in adata.var_names if g in genes_sp]]
    
    adata.write_h5ad(file_out)
