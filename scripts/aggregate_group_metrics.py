import pandas as pd
import numpy as np
import anndata as ad
import txsim.txsim as tx
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
    rand_matrix_list = []
    ann_matrix_list = []
    
    #Read in all the individual matrices
    for replicate in replicate_folders:
        matrix_file_name = os.path.join(replicate, f"rand_matrix.csv") # TODO FIX TO INCLUDE OTHER THINGS
        rand_matrix_list.append( pd.read_csv(matrix_file_name, index_col = 0 ))

        matrix_file_name = os.path.join(replicate, f"annotation_matrix.csv")
        ann_matrix_list.append( pd.read_csv(matrix_file_name, index_col = 0 ))

    [mean_rand_matrix, std_rand_matrix] = tx.metrics.aggregate_rand_index(rand_matrix_list)
    [mean_ann_matrix, std_ann_matrix] = tx.metrics.aggregate_rand_index(ann_matrix_list) #the rand aggregate also works for this so...

    #Save the mean and standard deviation matrices
    mean_rand_matrix.to_csv(f'{data}/aggregated/mean_rand_matrix.csv')
    std_rand_matrix.to_csv(f'{data}/aggregated/std_rand_matrix.csv')
    
    mean_ann_matrix.to_csv(f'{data}/aggregated/mean_annotation_matrix.csv')
    std_ann_matrix.to_csv(f'{data}/aggregated/std_annotation_matrix.csv')

    metric_list = pd.DataFrame()
    metric_list.to_csv(os.path.join(data, f"aggregated/aggregated_group_metrics.csv"))
    