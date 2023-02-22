import pandas as pd
import numpy as np
import anndata as ad
import txsim as tx
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate metrics across different runs')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Output data directory- should also count matrix')
    parser.add_argument('-f', '--files', required=True, type=str,
        help='List of all files or "all" for all runs in dataset')
    
    args = parser.parse_args()
    
    data = args.data
    files = args.files

    assignments = pd.DataFrame()

    if files == 'all':
        count_files = list( filter(lambda file: 'counts_' in file, os.listdir(data)))
    else:
        count_files = list(eval(files)) #TODO fix every instance of eval -> not good style apparently

    for count_matrix in count_files:
        path = os.path.join(data, count_matrix)
        adata = ad.read_h5ad(path)
        row_name = count_matrix.replace('counts_', '').replace('.h5ad','')
        assignments[row_name] = adata.uns['spots'].fillna(0).replace({-1:0})['cell']

    df = tx.metrics.calc_rand_index(assignments)
    df.to_csv(f'{data}/rand_matrix.csv')

    metric_list = pd.DataFrame()
    metric_list.to_csv(f'{data}/group_metrics.csv')

    