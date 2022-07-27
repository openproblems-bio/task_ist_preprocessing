#!/usr/bin/env python

import scanpy as sc
import txsim as tx
import argparse

#INPUT: (Spatial) counts.h5ad, (scRNAseq) counts.h5ad
#OUTPUT: metrics.txt
#From config
#sc_data = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/heart_sc.h5ad" 
#segmentation_method = 'imagej'
#assignment_method = 'pciSeq'
#normalize_by = 'area'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also count matrix')
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

    #Generate metrics
    coex_all = tx.metrics.coexpression_similarity(stdata, scdata)
    coex_thresh = tx.metrics.coexpression_similarity(stdata, scdata, thresh=0.5)

    f = open(f'{data}/metrics_{methods}.txt', 'w')
    f.truncate(0)
    f.writelines(f"Coexpresion Difference: {coex_all}\n")
    f.writelines(f"Coexpresion Difference, with threshold: {coex_thresh}")
    f.close()
