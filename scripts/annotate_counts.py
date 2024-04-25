#!/usr/bin/env python

import yaml
from pathlib import Path
import txsim as tx
import scanpy as sc
import pandas as pd
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate count matrix for spatial data')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain assignments_.csv')
    parser.add_argument('-s', '--singlecell', required=True, type=str,
	    help='Path to the single cell anndata')
    parser.add_argument('-c', '--counts', required=True, type=str, 
        help='Method list after counts_')
    parser.add_argument('-n', '--annotation_csv', required=True, type=str, 
        help='CSV file with cell type annotations')
    parser.add_argument('-a', '--annotate', required=True, type=str,
        help='Method to annotate celltypes') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Optional dictionary (as string) of method-specific parameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters') 

    args = parser.parse_args()

    counts_method = args.counts
    data = args.data
    ct_method = args.annotate
    id_code = args.id_code
    file_sc = args.singlecell

    # Get hyperparameters
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults[ct_method]
        gparams_defaults = defaults["annotation_params"] 
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    groupparams = eval(args.groupparams)
    groupparams.update({k:v for k,v in gparams_defaults.items() if k not in groupparams})
    groupparams = {k:(v if v != "None" else None) for k,v in groupparams.items()}
    
    # Variables for specific params
    per_gene_correction = groupparams['per_gene_correction']
    gene_corr_layer = groupparams['per_gene_layer']

    # Read in the single-cell data
    adata_sc = sc.read(file_sc)
    adata = sc.read(f'{data}/normcounts_{counts_method}.h5ad')
    
    # Read in the cell type annotations
    df_cts = pd.read_csv(args.annotation_csv, index_col=0)

    # Apply filters (TODO: keep/remove? extend to other methods?)
    if ct_method in ['majority', 'ssam', 'pciseqct']: 
        df_cts = df_cts.loc[df_cts["score"] < hyperparams["ct_threshold"], "celltype"] = "None_sp"

    # Transfer cell type annotations to spatial data
    adata.obs['celltype'] = df_cts.loc[adata.obs['cell_id'], "celltype"]

    # Do per-gene correction if active
    if per_gene_correction:
        tx.preprocessing.gene_efficiency_correction(adata, adata_sc, gene_corr_layer)
        if gene_corr_layer != 'lognorm':
            adata.layers['lognorm'] = adata.layers['norm']
            sc.pp.log1p(adata, layer='lognorm')

    #Save AnnData object
    adata.write_h5ad(f"{data}/counts_{counts_method}_{ct_method}-{id_code}.h5ad")
