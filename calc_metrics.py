#!/usr/bin/env python

import scanpy as sc
import txsim as tx

#INPUT: (Spatial) counts.h5ad, (scRNAseq) counts.h5ad
#OUTPUT: metrics.txt
#From config
sc_data = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/heart_sc.h5ad" 
segmentation_method = 'imagej'
assignment_method = 'pciSeq'
normalize_by = 'area'

#Read count matrices
scdata = sc.read(sc_data)
stdata = sc.read(f'data/counts_{segmentation_method}_{assignment_method}_{normalize_by}.h5ad')

#Generate metrics
coex_all = tx.metrics.coexpression_similarity(stdata, scdata)
coex_thresh = tx.metrics.coexpression_similarity(stdata, scdata, thresh=0.5)

f = open(f'data/metrics_{segmentation_method}_{assignment_method}_{normalize_by}', 'w')
f.truncate(0)
f.writelines(f"Coexpresion Difference: {coex_all}\n")
f.writelines(f"Coexpresion Difference, with threshold: {coex_thresh}")
f.close()

