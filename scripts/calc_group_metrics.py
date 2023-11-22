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
        help='Output data directory')
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    args = parser.parse_args()
    
    data = args.data
    id_code = args.id_code

    # join with parent since data include replicate
    df = pd.read_csv(os.path.join(os.path.dirname(data), 'group_metric_chunks.csv'), index_col = False)
    # read relevant pairs based on chunk id
    chunk_df = df.loc[df['chunk_id']==int(id_code),['run1','run2']]
    #print(chunk_df)
    index_list = list(chunk_df.itertuples(index=False))

    metric_list = pd.DataFrame(
        columns = ['rand_index', 'annotation_similarity'], 
        index = pd.MultiIndex.from_tuples(index_list, names=('run1','run2')))

    #Calculate each metric
    #TODO make more efficient by not reloading the same one every time
    for pair in index_list:
        adata1 = ad.read_h5ad(os.path.join(data, 'counts_' + pair[0] + '.h5ad'))
        adata2 = ad.read_h5ad(os.path.join(data, 'counts_' + pair[1] + '.h5ad'))
        
        ann_sim = tx.metrics.calc_annotation_similarity(
            adata1, adata2)

        #TODO always change to adjusted rand index?
        rand_idx = sklearn.metrics.adjusted_rand_score(
            adata1.uns['spots']['cell'].fillna(0).replace({-1:0}), 
            adata2.uns['spots']['cell'].fillna(0).replace({-1:0}))
        
        metric_list['rand_index'][pair] = rand_idx
        metric_list['annotation_similarity'][pair] = ann_sim
    
    metric_list.to_csv(f'{data}/group_metrics-{id_code}.csv')

    