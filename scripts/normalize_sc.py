#!/usr/bin/env python

import txsim as tx
import scanpy as sc
import os.path
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-sc', '--singlecell', required=True, type=str, 
        help='Non-normalized adata')
    parser.add_argument('-o', '--output', default='total', type=str,
        help='Adata output') 
    
    args = parser.parse_args()

    file_in = args.singlecell
    file_out = args.output

    adata = sc.read(file_in)
    adata = tx.preprocessing.normalize_sc(adata)
	
    adata.write_h5ad(file_out)
