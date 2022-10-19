#!/usr/bin/env python

import argparse
import txsim as tx
import pandas as pd
import skimage.io

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
 
    spots = tx.preprocessing.basic_assign(
        pd.read_csv(molecules),
        skimage.io.imread(f'{data}/segments_{segmentation_method}.tif')
    )

    #Save in correct format
    spots.to_csv(f'{data}/assignments_{segmentation_method}_basic-{id_code}.csv', index = False)