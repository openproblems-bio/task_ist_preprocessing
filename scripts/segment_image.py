#!/usr/bin/env python

import txsim as tx
import squidpy as sq
from PIL import Image
from scipy import ndimage
import skimage.io
import skimage.measure
import skimage.segmentation
import numpy as np
import argparse
import os

#TODO: change cellpose and watershed to 2 separate methods 

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
    parser.add_argument('-e', '--expand', default=0, type=int,
        help='Amount to expand each segment by- can be used to approximate cell boundary') 
    
    args = parser.parse_args()

    image_file = args.input
    output = args.output
    binary = args.binary
    segmentation_method = args.segment
    expand_nuclear_area = args.expand
    id_code = args.id_code
    
    #Create output folder if needed
    if not os.path.exists(output):
        os.mkdir(output)

    #If unsegmented, segment image
    if(not binary):
        img = sq.im.ImageContainer(image_file)
        tx.preprocessing.segment_nuclei(img, layer = 'image', method=segmentation_method)
        img_arr = img[f'segmented_{segmentation_method}'].to_numpy()[:,:,0,0]
    
    #If already segmented, label
    else:
        img_arr = skimage.io.imread(image_file)
        img_arr = skimage.measure.label(img_arr, connectivity=1)

    #Expand nuclear area to reflect whole cell area
    if(expand_nuclear_area != None):
        img_arr = skimage.segmentation.expand_labels(img_arr, distance=expand_nuclear_area)

    #Save as .tif file
    skimage.io.imsave(f'{output}/segments_{segmentation_method}-{id_code}.tif', img_arr)

    #Calculate and save areas
    (unique, counts) = np.unique(img_arr, return_counts=True)
    areas = np.asarray((unique, counts)).T
    np.savetxt(f'{output}/areas_{segmentation_method}-{id_code}.csv', areas, delimiter=",")
