#!/usr/bin/env python

import argparse
import pciSeq
import scipy.io
from scipy.sparse import coo_matrix
import pandas as pd
import scanpy as sc

#INPUT: molecules.csv, singlecell.h5ad labels.mat
#OUTPUT: assignments.csv
#From config:
#segmentation_method = 'imagej'
#molecules = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/spots_PCW4.5_1.csv"
#sc_data = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/heart_sc.h5ad" 
#opts = {'exclude_genes': ['TCIM']}

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-s', '--segment', required=True, type=str,
        help='Segmentation method used for image')
    parser.add_argument('-sc', '--singlecell', required=True, type=str,
        help='Single cell h5ad count matrix with celltype in anndata.obs[\'celltype\']') 
    parser.add_argument('-o', '--opts', default=None, type=dict,
        help='Option dictionary for pciSeq') 
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    segmentation_method = args.segment
    sc_data = args.singlecell
    opts = args.opts

    #Read and format molecules, single cell data, and labels
    spots = pd.read_csv(molecules)
    spots.columns = ['Gene', 'x', 'y']

    adata = sc.read_h5ad(sc_data)
    scdata = adata.X
    scdata  = pd.DataFrame(scdata.transpose())
    scdata.columns = adata.obs['celltype']
    scdata.index = adata.var_names

    label = scipy.io.loadmat(f'{data}/label_{segmentation_method}')['label']
    coo = coo_matrix(label)

    #Run through pciSeq
    pciSeq.attach_to_log()
    if(opts != None):
        cellData, geneData = pciSeq.fit(spots, coo, scdata, opts)
    else:
        cellData, geneData = pciSeq.fit(spots, coo, scdata)   

    #Save in correct format
    assignments = geneData[ ["Gene", "x", "y", "neighbour"] ]
    assignments.columns = ["gene", "x", "y", "cell"]
    assignments.to_csv(f'{data}/assignments_{segmentation_method}_pciseq.csv', index = False)
