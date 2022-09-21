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
    
    args = parser.parse_args()
    
    methods = args.methods
    data = args.data

    #Read count matrices
    stdata = sc.read(f'{data}/counts_{methods}.h5ad')

    #Gather information on spatial data
    df = tx.quality_metrics.all_quality_metrics(stdata)

    df.to_csv(f'{data}/quality_metrics_{methods}.csv')
    
