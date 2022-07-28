#!/usr/bin/env python

import argparse
import pciSeq
import skimage.io
from scipy.sparse import coo_matrix
import pandas as pd
import scanpy as sc
import numpy as np

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
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    segmentation_method = args.segment
    sc_data = args.singlecell
    opts = eval(args.hyperparams)['opts']
    id_code = args.id_code
 
    #Read and format molecules, single cell data, and labels
    spots = pd.read_csv(molecules)
    spots.columns = ['Gene', 'x', 'y']

    adata = sc.read_h5ad(sc_data)
    scdata = adata.X
    scdata  = pd.DataFrame(scdata.transpose())
    scdata.columns = adata.obs['celltype']
    scdata.index = adata.var_names

    seg = skimage.io.imread(f'{data}/segments_{segmentation_method}.tif')
    coo = coo_matrix(seg)

    #Add safety feature for genes that aren't included
    #Run through pciSeq
    pciSeq.attach_to_log()
    if(opts != None):
        cellData, geneData = pciSeq.fit(spots, coo, scdata, opts)
    else:
        cellData, geneData = pciSeq.fit(spots, coo, scdata)   

    #Save in correct format
    assignments = geneData[ ["Gene", "x", "y", "neighbour"] ]
    assignments.columns = ["gene", "x", "y", "cell"]

    #Change the cell names to match the segmentation
    assignments['cell'] = np.unique(seg)[assignments['cell']]
    assignments.to_csv(f'{data}/assignments_{segmentation_method}_pciseq-{id_code}.csv', index = False)
