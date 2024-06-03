#!/usr/bin/env python

from pathlib import Path
import argparse
from collections import OrderedDict
import numpy as np
import pandas as pd
import anndata as ad
import txsim as tx
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
    
    # Load image and assign cell ids such that there are no gaps 
    # (NOTE, TODO: Arbitrary ids lead to an error in tx.preprocessing.run_pciSeq. The most elegant solution would be to 
    # keep the ids. But in general thr requirements wrt cell_ids in the assignment outputs need to be checked. Also, 
    # pciseq needs to be revisited in the txsim package under consideration of a stable release tag. It seems that 
    # there were quite some differences wrt the methods outputs.)
    image = skimage.io.imread(f'{data}/segments_{segmentation_method}.ome.tif')
    background_exists = 0 in image
    image = np.unique(image, return_inverse=True)[1].reshape(image.shape)
    if not background_exists:
        image += 1
        image[0,0] = 0 # introduce 1 pixel as pseudo background
    
    #Read data and run pciSeq
    assignments, cell_types = tx.preprocessing.run_pciSeq(
        pd.read_csv(molecules),
        image,
        ad.read(sc_data),
        'celltype',
        opts
    )

    #Save to csv
    cell_types.to_csv(f'{data}/celltypes_{segmentation_method}_pciseq-{id_code}.csv')
    assignments.to_csv(f'{data}/assignments_{segmentation_method}_pciseq-{id_code}.csv', index = False)
