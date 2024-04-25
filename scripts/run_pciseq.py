#!/usr/bin/env python

from pathlib import Path
import argparse
from collections import OrderedDict
import txsim as tx
import pandas as pd
import anndata as ad
import skimage.io
import yaml


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
    

    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        hparams_defaults = yaml.safe_load(file)["pciseq"]
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    segmentation_method = args.segment
    sc_data = args.singlecell
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    opts = dict(hyperparams.get('opts')) if (hyperparams is not None and hyperparams.get('opts') is not None) else None
    id_code = args.id_code
    
    #Read data and run pciSeq
    assignments, cell_types = tx.preprocessing.run_pciSeq(
        pd.read_csv(molecules),
        skimage.io.imread(f'{data}/segments_{segmentation_method}.ome.tif'),
        ad.read(sc_data),
        'celltype',
        opts
    )

    #Save to csv
    cell_types.to_csv(f'{data}/celltypes_{segmentation_method}_pciseq-{id_code}.csv')
    assignments.to_csv(f'{data}/assignments_{segmentation_method}_pciseq-{id_code}.csv', index = False)
