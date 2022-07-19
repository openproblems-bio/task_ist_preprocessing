#!/usr/bin/env python
# coding: utf-8

import txsim as tx
import squidpy as sq
import scipy.io
from scipy import ndimage
import skimage.io
import skimage.measure
import skimage.segmentation
import numpy as np

#INPUT: DAPI.tif, binary.tif
#OUTPUT: label.mat, areas.csv
#From config:
image_file = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/nuclei_PCW4.5_1_watershed.tif"
binary = True #Note: assumes the binary image is segmented (e.g. watershed via ImageJ)
segmentation_method = 'imagej'
expand_nuclear_area = 10

#If unsegmented, segment image
if(not binary):
    img = sq.im.ImageContainer(image_file)
    tx.preprocessing.segment_nuclei(img, layer = 'image', method=segmentation_method)
    img_arr = img["segmented_{}".format(segmentation_method)].to_numpy()[:,:,0,0]

#If already segmented, label
else:
    img_arr = skimage.io.imread(image_file)
    img_arr = skimage.measure.label(img_arr, connectivity=1)

#Expand nuclear area to reflect whole cell area
if(expand_nuclear_area != None):
    img_arr = skimage.segmentation.expand_labels(img_arr, distance=expand_nuclear_area)

#Save as .mat file
scipy.io.savemat("data/label_{}.mat".format(segmentation_method), {'label':img_arr})

#Calculate and save areas
(unique, counts) = np.unique(img_arr, return_counts=True)
areas = np.asarray((unique, counts)).T
np.savetxt('data/areas_{}.csv'.format(segmentation_method), areas, delimiter=",")

