import pandas as pd
import numpy as np
import anndata as ad
import txsim as tx
import argparse
import os
import itertools
import sklearn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate metrics across different runs')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Output data directory- should also count matrix')
    parser.add_argument('-f', '--files', required=True, type=str,
        help='List of all files or "all" for all runs in dataset')
    
    args = parser.parse_args()
    
    data = args.data
    files = args.files

    if files == 'all':
        count_files = list( filter(lambda file: 'counts_' in file and 'norm' not in file, os.listdir(data)))
        normcount_files = list( filter(lambda file: 'normcounts_' in file, os.listdir(data)))
    else:
        count_files = list( filter(lambda file: 'counts_' in file and 'norm' not in file, eval(files)))
        normcount_files = list( filter(lambda file: 'normcounts_' in file, eval(files))) #TODO fix every instance of eval -> not good style apparently
    
    index_list = []
    for pair in itertools.combinations(count_files,2):
        name1 = pair[0].replace('counts_', '').replace('.h5ad','')
        name2 = pair[1].replace('counts_', '').replace('.h5ad','')
        index_list.append(frozenset([name1,name2]))
    index_list = [tuple(name) for name in list(set(index_list))]

    metric_list = pd.DataFrame(
        columns = ['rand_index', 'annotation_similarity'], 
        index = pd.MultiIndex.from_tuples(index_list))

    #Calculate each metric
    for pair in index_list:
        adata1 = ad.read_h5ad(os.path.join(data, 'counts_' + pair[0] + '.h5ad'))
        adata2 = ad.read_h5ad(os.path.join(data, 'counts_' + pair[1] + '.h5ad'))
        
        ann_sim = tx.metrics.calc_annotation_similarity(
            adata1, adata2)
    
        rand_idx = sklearn.metrics.rand_score(
            adata1.uns['spots']['cell'].fillna(0).replace({-1:0}), 
            adata2.uns['spots']['cell'].fillna(0).replace({-1:0}))
        
        metric_list['rand_index'][pair] = rand_idx
        metric_list['annotation_similarity'][pair] = ann_sim
    
    metric_list.to_csv(f'{data}/group_metrics.csv')

    