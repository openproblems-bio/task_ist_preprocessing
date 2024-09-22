import txsim as tx
#from pathlib import Path
#import tifffile
#import squidpy as sq
#from scipy import ndimage
#import skimage.io
#import skimage.measure
#import skimage.segmentation
import numpy as np
import txsim as tx 
#import argparse
import os
import yaml
import spatialdata as sd
import anndata as ad
import shutil
import numpy as np
from spatialdata.models import Labels2DModel
import xarray as xr
import datatree as dt


def convert_to_lower_dtype(arr):
    max_val = arr.max()
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64

    return arr.astype(new_dtype)

## VIASH START
par = {
  "input": "../task_ist_preprocessing/resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr",
  "output": "segmentation.zarr",
  'hyperparams': [] ##############################how to take from viash
}

## VIASH END

sdata = sd.read_zarr(par["input"])
image = sdata['rep1_image']['scale0'].image.compute().to_numpy()
transformation = image.transform.copy()
transformation['global'] = transformation.pop('rep1_global')
image = convert_to_lower_dtype(image)

sd_output = sd.SpatialData()
scales = [2000, 1000, 500, 250, 125]
downsampled_arrays = {}
for idx, scale in enumerate(scales):
    img_arr = tx.preprocessing.segment_cellpose(image, hyperparams)  
    data_array = xr.DataArray(img_arr, name=f'segmentation_scale{scale}', dims=('y', 'x'))
    parsed_data = Labels2DModel.parse(data_array, transformations=transformation)
    downsampled_arrays[f'scale{idx}'] = parsed_data
tree = dt.DataTree()
for scale_key, parsed_data in downsampled_arrays.items():
    tree[scale_key] = dt.DataTree(parsed_data)

sd_output.labels['segmentation'] = tree


print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sd_output.write(par["output"])

_____________________________________________






args -> put into par dict






(with tiffle) line   - we don't write in tiff file anymore, we convert the output to the spatialdata: coodinates with transformation, we have a crop image (I want to keep the transformation of the crop)
don;t do save area (at the bottom) (edited) 
_______________________________________________
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







___________________________________________________________________________
import spatialdata as sd
import anndata as ad
import os
import shutil

## VIASH START
par = {
  "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
  "labels_key": "cell_labels",
  "output": "segmentation.zarr",
}
meta = {
  "name": "segmentation"
}
## VIASH END

print("Reading input files", flush=True)
sdata = sd.read_zarr(par["input"])

assert par["labels_key"] in sdata.labels, f"Key '{par['labels_key']}' not found in input data."

print(f"Copy segmentation from '{par['labels_key']}'", flush=True)
sdata_segmentation_only = sd.SpatialData(
  labels={
    "segmentation": sdata[par["labels_key"]]
  },
  tables={
    "table": ad.AnnData(
      obs=sdata.tables["table"].obs[["cell_id", "region"]],
      var=sdata.tables["table"].var[[]]
    )
  }
)

print("Writing output", flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_segmentation_only.write(par["output"])
