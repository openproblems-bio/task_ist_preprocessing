import argparse
import os
from pathlib import Path
from deepcell.utils.io_utils import get_image
import numpy as np
from deepcell.applications import Mesmer
import tifffile
import yaml
import skimage.io

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


# Prepare input https://github.com/vanvalenlab/deepcell-applications/blob/master/deepcell_applications/prepare.py
# load the input files into numpy arrays
def load_image(path, channel=0, ndim=3):
    """Load an image file as a single-channel numpy array.
    Args:
        path (str): Filepath to the image file to load.
        channel (list): Loads the given channel if available.
            If channel is list of length > 1, each channel
            will be summed.
        ndim (int): The expected rank of the returned tensor.
    Returns:
        numpy.array: The image channel loaded as an array.
    """
    if not path:
        raise IOError('Invalid path: %s' % path)

    img = get_image(path)

    channel = channel if isinstance(channel, (list, tuple)) else [channel]

    # getting a little tricky, which axis is channel axis?
    if img.ndim == ndim:
        # file includes channels, find the channel axis
        # assuming the channels axis is the smallest dimension
        axis = img.shape.index(min(img.shape))
        if max(channel) >= img.shape[axis]:
            raise ValueError('Channel {} was passed but channel axis is '
                             'only size {}'.format(
                                 max(channel), img.shape[axis]))

        # slice out only the required channel
        slc = [slice(None)] * len(img.shape)
        # use integer to select only the relevant channels
        slc[axis] = channel
        img = img[tuple(slc)]
        # sum on the channel axis
        img = img.sum(axis=axis)

    # expand the (proper) channel axis
    img = np.expand_dims(img, axis=-1)

    if not img.ndim == ndim:
        raise ValueError('Expected image with ndim = {} but found ndim={} '
                         'and shape={}'.format(ndim, img.ndim, img.shape))

    return img



if __name__ == '__main__':
# Parse CLI taken from https://github.com/vanvalenlab/deepcell-applications/blob/master/deepcell_applications/argparse.py
    root_dir = os.path.dirname(os.path.abspath(__file__))
    
    class WritableDirectoryAction(argparse.Action):
        # Check that the provided output directory exists, and is writable
        # From: https://gist.github.com/Tatsh/cc7e7217ae21745eb181
        def __call__(self, parser, namespace, values, option_string=None):
            prospective_dir = values
            if not os.path.isdir(prospective_dir):
                raise argparse.ArgumentTypeError(
                    '{} is not a valid directory'.format(prospective_dir, ))
            if os.access(prospective_dir, os.W_OK | os.X_OK):
                setattr(namespace, self.dest, os.path.realpath(prospective_dir))
                return
            raise argparse.ArgumentTypeError(
                '{} is not a writable directory'.format(prospective_dir))

    def existing_file(x):
        if x is not None and not os.path.exists(x):
            raise argparse.ArgumentTypeError('{} does not exist.'.format(x))
        return x
    
    parser = argparse.ArgumentParser(description='Segment nuclear and optionally membrane stained cells with Mesmer')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output directory')
    #parser.add_argument('--output-directory', '-o',
    #                    default=os.path.join(root_dir, 'output'),
    #                    action=WritableDirectoryAction,
    #                    help='Directory where application outputs are saved.',dest='output_directory')
    #parser.add_argument('--output-name', '-f',
    #                    default='mask.tif',
    #                    help='Name of output file.')
    #parser.add_argument('--nuclear-image', '-i', required=True,
    #                    type=existing_file, dest='nuclear_path',
    #                    help=('REQUIRED: Path to 2D single channel TIF file.'))    
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-g', '--groupparams', default=None, type=str,
        help='Optional dictionary (as string) of group parameters') 
    
    # Get arguments
    args = parser.parse_args()
    
    image_file = args.input
    output = args.output
    id_code = args.id_code
    
    # Get hyperparameters
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults["mesmer"]
        gparams_defaults = defaults["segmentation_params"]
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    groupparams = eval(args.groupparams)
    groupparams.update({k:v for k,v in gparams_defaults.items() if k not in groupparams})
    groupparams = {k:(v if v != "None" else None) for k,v in groupparams.items()}
    expand_nuclear_area = groupparams.get('expand') #If None, it will not expand after segmenting    
    
    # Create specified output directory if it does not exist
    if not os.path.exists(output):
        os.makedirs(output)
    
    # Takes the user-supplied command line arguments and runs the specified application https://github.com/vanvalenlab/deepcell-applications/blob/master/deepcell_applications/app_runners.py
    
    # # Check that the output path does not exist already
    # if os.path.exists(outfile):
    #     raise IOError(f'{outfile} already exists!')
        
    # Load in the nuclear image
    nuclear_img = load_image(image_file, channel=0, ndim=3) #args.nuclear_channel
    
    # Optionally, load in the membrane image
    # if args.membrane_path:
    #         membrane_img = load_image(args.membrane_path, channel=args.membrane_channel, ndim=3)
    # else:
    membrane_img = np.zeros(nuclear_img.shape, dtype=nuclear_img.dtype)

    # join the inputs in the correct order
    img = np.concatenate([nuclear_img, membrane_img], axis=-1)
    

    # make sure the input image is compatible with the app
    # Grabbed from https://deepcell.readthedocs.io/en/master/API/deepcell.applications.html#mesmer to simplify
    model_image_shape = (128, 128, 1)
    rank = len(model_image_shape)
    errtext = ("Invalid image shape")

    if len(img.shape) != len(model_image_shape):
        raise ValueError(errtext)

    # Applications expect a batch dimension
    img = np.expand_dims(img, axis=0)
    
    # Create the application
    print("Creating Mesmer app")
    app = Mesmer()

    # run the prediction
    batch_size = hyperparams.get('batch_size')
    image_mpp = hyperparams.get('image_mpp')
    compartment = hyperparams.get('compartment')
    img_arr = app.predict(img, batch_size=batch_size,image_mpp=image_mpp,compartment=compartment)
    print("Prediction compelete")

    # Optionally squeeze the output
    if hyperparams.get('squeeze') is None or not hyperparams.get('squeeze'):
        img_arr = np.squeeze(img_arr)
        
    # Expand
    if expand_nuclear_area is not None and expand_nuclear_area != 0:
        img_arr = skimage.segmentation.expand_labels(img_arr, distance=expand_nuclear_area)

    # Save as .ome.tif file
    img_arr = convert_to_lower_dtype(img_arr)
    with tifffile.TiffWriter(f'{output}/segments_mesmer-{id_code}.ome.tif', bigtiff=True) as tif:
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
    tifffile.imwrite(f'{output}/segments_mesmer-{id_code}.tif', img_arr)

    #Calculate and save areas
    (unique, counts) = np.unique(img_arr, return_counts=True)
    areas = np.asarray((unique, counts)).T
    np.savetxt(f'{output}/areas_mesmer-{id_code}.csv', areas, delimiter=",")