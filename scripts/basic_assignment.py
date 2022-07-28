#!/usr/bin/env python

import argparse
import pciSeq
import skimage.io
import pandas as pd
import numpy as np
import scanpy as sc

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-s', '--segment', required=True, type=str,
        help='Segmentation method used for image')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    segmentation_method = args.segment
    id_code = args.id_code
    hyperparams = args.hyperparams
 
    #Read and format molecules
    spots = pd.read_csv(molecules)
    spots.columns = ['gene', 'x', 'y']
    seg = skimage.io.imread(f'{data}/segments_{segmentation_method}.tif')

    #Assign molecules based value in segmentation array
    spots['cell'] = seg[spots.y.to_numpy(dtype=np.int64), spots.x.to_numpy(dtype=np.int64)]

    #Save in correct format
    spots.to_csv(f'{data}/assignments_{segmentation_method}_basic-{id_code}.csv', index = False)