import pandas as pd
import numpy as np
import anndata as ad
import txsim as tx
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate metrics across different runs')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Output data directory- ')
    
    
    args = parser.parse_args()
    
    data = args.data

    if not os.path.exists( os.path.join(data, "aggregated") ):
        os.makedirs( os.path.join(data, "aggregated") )

    replicate_folders = [f.path for f in os.scandir(data) if f.is_dir() and 'replicate' in f.path]
    matrix_list = []
    
    for replicate in replicate_folders:
        matrix_file_name = os.path.join(replicate, f"rand_matrix.csv") # TODO FIX TO INCLUDE OTHER THINGS
        matrix_list.append( pd.read_csv(matrix_file_name, index_col = 0 ))

    [mean_matrix, std_matrix] = tx.metrics.aggregate_rand_index(matrix_list)

    mean_matrix.to_csv(f'{data}/aggregated/mean_rand_matrix.csv')
    std_matrix.to_csv(f'{data}/aggregated/std_rand_matrix.csv')
    
    metric_list = pd.DataFrame()
    metric_list.to_csv(os.path.join(data, f"aggregated/aggregated_group_metrics.csv"))
    