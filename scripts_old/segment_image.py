#!/usr/bin/env python

import txsim as tx
from pathlib import Path
import tifffile
import squidpy as sq
from scipy import ndimage
import skimage.io
import skimage.measure
import skimage.segmentation
import numpy as np
import argparse
import os
import yaml

def convert_to_lower_dtype(arr):
    # Find the maximum value in the array
    max_val = arr.max()

    # Determine the smallest unsigned integer dtype that can hold the max value
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64

    # Convert the array to the determined dtype
    return arr.astype(new_dtype)


if __name__ == '__main__':

    #Parse arguments
    parser = argparse.ArgumentParser(description='Segment DAPI images')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output directory')
    parser.add_argument('-b', '--binary', action='store_true',
        help='If the input image is a segmented, binary image (e.g. watershed via ImageJ)')
    parser.add_argument('-s', '--segment', required=True, type=str,
        help='Segmentation method to be used')
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters') 
    
    args = parser.parse_args()

    image_file = args.input
    output = args.output
    binary = args.binary
    segmentation_method = args.segment
    id_code = args.id_code
    
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults[segmentation_method]
        gparams_defaults = defaults["segmentation_params"]
    
    print("################## ", args.groupparams)
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    groupparams = eval(args.groupparams)
    groupparams.update({k:v for k,v in gparams_defaults.items() if k not in groupparams})
    groupparams = {k:(v if v != "None" else None) for k,v in groupparams.items()}
    expand_nuclear_area = groupparams.get('expand') #If None, it will not expand after segmenting
    
    #Create output folder if needed
    if not os.path.exists(output):
        os.makedirs(output)

    #If unsegmented, segment image
    if(not binary):
        if(segmentation_method == 'binning'):
            img = tifffile.imread(image_file)
            #if hyperparams is None or hyperparams['bin_size'] is None:
            #    img_arr = tx.preprocessing.segment_binning(img, 20)
            #else:
            img_arr = tx.preprocessing.segment_binning(img, hyperparams['bin_size'])
        elif(segmentation_method == 'stardist'):
            img = tifffile.imread(image_file)
            img_arr = tx.preprocessing.segment_stardist(img, hyperparams)
        elif(segmentation_method=='cellpose'):
            img = tifffile.imread(image_file)
            print("In cellpose!")
            img_arr = tx.preprocessing.segment_cellpose(img, hyperparams)

        elif(segmentation_method=='watershed'):
            print("In watershed!")
            img = tifffile.imread(image_file)
            img_arr = tx.preprocessing.segment_watershed(img, hyperparams)
            
        else: #TODO: Do we still need this?
            img = sq.im.ImageContainer(image_file)
            if hyperparams is not None:
                tx.preprocessing.segment_nuclei(img, layer = 'image', method=segmentation_method, **hyperparams)
            else:
                tx.preprocessing.segment_nuclei(img, layer = 'image', method=segmentation_method)
            img_arr = img[f'segmented_{segmentation_method}'].to_numpy()[:,:,0,0]
    
    else: #If already segmented, label
        img_arr = skimage.io.imread(image_file)
        img_arr = skimage.measure.label(img_arr, connectivity=1)
    print("Segmentation Done!")
    #Expand nuclear area to reflect whole cell area
    if expand_nuclear_area is not None and expand_nuclear_area != 0:
        img_arr = skimage.segmentation.expand_labels(img_arr, distance=expand_nuclear_area)
    
    img_arr = convert_to_lower_dtype(img_arr)
    #Save as .tif file
    #skimage.io.imsave(f'{output}/segments_{segmentation_method}-{id_code}.tif', img_arr)
    with tifffile.TiffWriter(f'{output}/segments_{segmentation_method}-{id_code}.ome.tif', bigtiff=True) as tif:
        metadata={
            'PhysicalSizeX': 1,
            'PhysicalSizeXUnit': 'um',
            'PhysicalSizeY': 1,
            'PhysicalSizeYUnit': 'um'
        }
        tif.write(
            img_arr,
            metadata=metadata
        )
    tifffile.imwrite(f'{output}/segments_{segmentation_method}-{id_code}.tif', img_arr)

    #Calculate and save areas
    (unique, counts) = np.unique(img_arr, return_counts=True)
    areas = np.asarray((unique, counts)).T
    np.savetxt(f'{output}/areas_{segmentation_method}-{id_code}.csv', areas, delimiter=",")
