#!/usr/bin/env python

import argparse
from collections import OrderedDict
import txsim as tx

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
    hyperparams = eval(args.hyperparams)
    opts = dict(hyperparams.get('opts')) if hyperparams is not None else None
    id_code = args.id_code
 
    assignments, cell_types = tx.preprocessing.run_pciSeq(
        molecules,
        f'{data}/segments_{segmentation_method}.tif',
        sc_data,
        'celltype',
        opts
    )

    #Save to csv
    cell_types.to_csv(f'{data}/celltypes_{segmentation_method}_pciSeq-{id_code}.csv', header = False)
    assignments.to_csv(f'{data}/assignments_{segmentation_method}_pciSeq-{id_code}.csv', index = False)
