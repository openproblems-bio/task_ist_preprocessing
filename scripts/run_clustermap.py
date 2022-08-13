#!/usr/bin/env python

import argparse
from collections import OrderedDict
import txsim as tx

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells using ClusterMap')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    image = args.image
    hyperparams = eval(args.hyperparams)
    id_code = args.id_code
 
    #TODO Add all of ClusterMap's many parameters
    
    assignments = tx.preprocessing.run_clustermap(
        molecules,
        image,
    )

    #Save to csv
    assignments.to_csv(f'{data}/assignments_clustermap-{id_code}.csv', index = False)
