#!/usr/bin/env python
import pandas as pd
import numpy as np
import scanpy as sc
import txsim as tx
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Output data directory- should also contain count matrix')
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
    df = tx.metrics.all_metrics(stdata, scdata)
    df.to_csv(f'{data}/metrics_{methods}.csv')

    if not os.path.exists(f'{data}/filtered'):
        os.makedirs(f'{data}/filtered')

    #Also with filtered cells
    df_filtered = tx.metrics.all_metrics(stdata[stdata.obs['passed_QC']], scdata)
    df_filtered.to_csv(f'{data}/filtered/metrics_{methods}.csv')
    