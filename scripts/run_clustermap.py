#!/usr/bin/env python

import argparse
from collections import OrderedDict
from ClusterMap.clustermap import ClusterMap
import tifffile
import numpy as np
import pandas as pd

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
    image = args.input
    hyperparams = eval(args.hyperparams)
    id_code = args.id_code

    #Make sure parameter dictionaries are not `None`
    if hyperparams is None: hyperparams = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    #gene_spots = load_gene_spot(gene_spot_path, folder_name, adjust_coords)
    #gene_spots.columns=['gene_name','spot_location_1','spot_location_2']
    ### convert gene_name to gene identity
    #genes=pd.DataFrame(gene_spots['gene_name'].unique())
    #a1=list(genes[0])
    #gene=list(map(lambda x: a1.index(x)+1, gene_spots['gene_name']))
    #gene_spots['gene']=gene
    #gene_spots['gene']=gene_spots['gene'].astype('int')
    
    xy_radius = hyperparams.get("model_xy_radius",40)
    z_radius = hyperparams.get("model_z_radius",0)
    fast_preprocess = hyperparams.get("model_fast_preprocess",False)
    window_size = hyperparams.get("model_window_size", 40)
    gauss_blur = hyperparams.get("model_gauss_blur", True)
    sigma = hyperparams.get("model_sigma", 1)
    trim =  hyperparams.get("model_trim", False)
    use_dapi =   hyperparams.get("model_use_dapi", True)
    
    window_size = hyperparams.get("preprocess_window_size", 1500)
    

    dapi_grid_interval = hyperparams.get("preprocess_dapi_grid_interval",5)
    contamination = hyperparams.get("preprocess_contamination", 0)
    pct_filter = hyperparams.get("preprocess_pct_filter", 0.05)
    LOF = hyperparams.get("preprocess_LOF", False)

    min_spot_per_cell = hyperparams.get("segmentation_min_spot_per_cell",5)
    
    
    cell_num_threshold = hyperparams.get("segmentation_cell_num_threshold", 0.01)
    dapi_grid_interval = hyperparams.get("segmentation_dapi_grid_interval", 5)
    add_dapi = hyperparams.get("segmentation_add_dapi", True)
    use_genedis = hyperparams.get("segmentation_use_genedis", True)
    
    #save gene annotation as genelist.csv
    #genes.to_csv(gene_spot_path+'_genelist.csv', header=None, index=False)
    #num_gene=np.max(gene_spots['gene'])

    #gene_list=np.arange(1,num_gene+1)
    #num_dims=len(img.shape)

    #Read and format input data
    dapi = tifffile.imread(image)
    num_dims=len(dapi.shape)
    spots = pd.read_csv(molecules)

    spots.rename(columns = {spots.columns[0]:'gene_name',
                        spots.columns[1]:'spot_location_1',
                        spots.columns[2]:'spot_location_2'
                        } , inplace = True)

    #Use gene id numbers instead of names
    genes, ids = np.unique(spots['gene_name'], return_inverse=True)
    spots['gene'] = ids+1
    spots = spots.astype({'spot_location_1':int, 'spot_location_2':int})
    gene_list=np.unique(ids)+1
    genes = pd.DataFrame(genes)

    #Create Model
    if use_dapi:
        if not trim:
            
            model = ClusterMap(spots=spots, dapi=dapi, gene_list=gene_list, num_dims=num_dims,
                xy_radius=xy_radius,z_radius=z_radius,fast_preprocess=fast_preprocess,gauss_blur=gauss_blur,
                sigma=sigma)

            #TODO: noise processing?
            ###preprocessing
            model.preprocess(dapi_grid_interval=dapi_grid_interval,pct_filter=pct_filter)
            ### segmentation
            model.min_spot_per_cell=min_spot_per_cell
            model.segmentation(cell_num_threshold=cell_num_threshold,
                                    dapi_grid_interval=dapi_grid_interval,add_dapi=add_dapi,use_genedis=use_genedis)
#print('Radius: ' + str(model.xy_radius))

    #The original spots file is modified by Clustermap to include assignments
    #Copy and return spots-to-cell assignment
    assignments = spots.copy()    
    assignments.rename(columns = {'gene_name':'Gene',
                        'spot_location_1':'x',
                        'spot_location_2':'y',
                        'clustermap':'cell',
                        } , inplace = True)
    assignments.cell = assignments.cell+1

    #Save to csv
    assignments.to_csv(f'{data}/assignments_clustermap-{id_code}.csv', index = False)
