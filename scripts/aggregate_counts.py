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
    parser.add_argument('-m', '--methods', required=True, type=str,
        help='Methods used to generate count matrices')
    
    args = parser.parse_args()
    
    data = args.data
    methods = args.methods

    if not os.path.exists( os.path.join(data, "aggregated") ):
        os.makedirs( os.path.join(data, "aggregated") )

    replicate_folders = [ f.path for f in os.scandir(data) if f.is_dir() and 'replicate' in f.path]
    adata_list = []
    rep_list = []
    
    for replicate in replicate_folders:
        count_file_name = os.path.join(replicate, f"counts_{methods}.h5ad")
        adata_list.append( ad.read(count_file_name) )
        rep_list.append( int(replicate.split("/")[-1].replace("replicate", "")) )


    adata = tx.preprocessing.aggregate_count_matrices(adata_list, rep_list)


    adata.write_h5ad(os.path.join(data, f"aggregated/counts_{methods}.h5ad"))
    