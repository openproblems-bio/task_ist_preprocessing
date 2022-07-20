#!/usr/bin/env python

import txsim as tx
import squidpy as sq
import scipy.io
from scipy import ndimage
import skimage.io
import skimage.measure
import skimage.segmentation
import numpy as np
import argparse

#INPUT: DAPI.tif, binary.tif
#OUTPUT: label.mat, areas.csv
#From config:
#image_file = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/nuclei_PCW4.5_1_watershed.tif"
#binary = True
#segmentation_method = 'imagej'
#expand_nuclear_area = 10


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Segment DAPI images')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output directory')
    parser.add_argument('-b', '--binary', default=False, type=bool,
        help='If the input image is a segmented, binary image (e.g. watershed via ImageJ)')
    parser.add_argument('-s', '--segment', default='watershed', type=str,
        help='Segmentation method to be used')
    parser.add_argument('-e', '--expand', default=0, type=int,
        help='Amount to expand each segment by- can be used to approximate cell boundary') 
    
    args = parser.parse_args()

    image_file = args.input
    output = args.output
    binary = args.binary
    segmentation_method = args.segment
    expand_nuclear_area = args.expand

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

    #Save as .mat file
    scipy.io.savemat(f'{output}/label_{segmentation_method}.mat', {'label':img_arr})

    #Calculate and save areas
    (unique, counts) = np.unique(img_arr, return_counts=True)
    areas = np.asarray((unique, counts)).T
    np.savetxt(f'{output}/areas_{segmentation_method}.csv', areas, delimiter=",")
