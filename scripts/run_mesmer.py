import argparse
import os
from deepcell.utils.io_utils import get_image
import numpy as np
from deepcell.applications import Mesmer
import tifffile
import skimage.io

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
                    '{} is not a valid directory'.format(
                        prospective_dir, ))
            if os.access(prospective_dir, os.W_OK | os.X_OK):
                setattr(namespace, self.dest, os.path.realpath(prospective_dir))
                return
            raise argparse.ArgumentTypeError(
                '{} is not a writable directory'.format(
                    prospective_dir))

    def existing_file(x):
        if x is not None and not os.path.exists(x):
            raise argparse.ArgumentTypeError('{} does not exist.'.format(x))
        return x
    
    parser = argparse.ArgumentParser(description='Segment nuclear and optionally membrane stained cells with Mesmer')
    
    parser.add_argument('--output-directory', '-o',
                        default=os.path.join(root_dir, 'output'),
                        action=WritableDirectoryAction,
                        help='Directory where application outputs are saved.',dest='output_directory')

    parser.add_argument('--output-name', '-f',
                        default='mask.tif',
                        help='Name of output file.')

    # parser.add_argument('--squeeze', action='store_true',
    #                     help='Squeeze the output tensor before saving.')

    parser.add_argument('--nuclear-image', '-i', required=True,
                        type=existing_file, dest='nuclear_path',
                        help=('REQUIRED: Path to 2D single channel TIF file.'))

    # parser.add_argument('--nuclear-channel', '-nc',
    #                     default=0, nargs='+', type=int,
    #                     help='Channel(s) to use of the nuclear image. '
    #                          'If more than one channel is passed, '
    #                          'all channels will be summed.')

    # parser.add_argument('--membrane-image', '-m',
    #                     type=existing_file, dest='membrane_path',
    #                     help=('Path to 2D single channel TIF file. '
    #                           'Optional. If not provided, membrane '
    #                           'channel input to network is blank.'))

    # parser.add_argument('--membrane-channel', '-mc',
    #                     default=0, nargs='+', type=int,
    #                     help='Channel(s) to use of the membrane image. '
    #                          'If more than one channel is passed, '
    #                          'all channels will be summed.')

    # parser.add_argument('--image-mpp', type=float, default=0.5,
    #                     help='Input image resolution in microns-per-pixel. '
    #                          'Default value of 0.5 corresponds to a 20x zoom.')

    # parser.add_argument('--batch-size', '-b', default=4, type=int,
    #                     help='Batch size for `model.predict`.')

    # parser.add_argument('--compartment', '-c', default='whole-cell',
    #                     choices=('nuclear', 'whole-cell', 'both'),
    #                     help='The cellular compartment to segment.')
    
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-e', '--expand', default="0", type=str,
        help='Amount to expand each segment by- can be used to approximate cell boundary') 
    
    args = parser.parse_args()
    hyperparams = eval(args.hyperparams) if args.hyperparams is not None else {}
    arg_dict = vars(args)
    
    # Create specified output directory if it does not exist
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    
    # Takes the user-supplied command line arguments and runs the specified application https://github.com/vanvalenlab/deepcell-applications/blob/master/deepcell_applications/app_runners.py
    
    
    
    # # Check that the output path does not exist already
    # if os.path.exists(outfile):
    #     raise IOError(f'{outfile} already exists!')
        
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
    
    # Load in the nuclear image
    nuclear_img = load_image(args.nuclear_path, channel=0, ndim=3) #args.nuclear_channel
    
    # Optionally, load in the membrane omage
    # if args.membrane_path:
    #         membrane_img = load_image(
    #         args.membrane_path,
    #         channel=args.membrane_channel,
    #         ndim=3)
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
    batch_size = hyperparams.get('batch_size') if hyperparams.get('batch_size') is not None else 1
    image_mpp = hyperparams.get('image_mpp') if hyperparams.get('image_mpp') is not None else 0.5
    compartment = hyperparams.get('compartment') if hyperparams.get('compartment') is not None else "nuclear"
    img_arr = app.predict(img, batch_size=batch_size,image_mpp=image_mpp,compartment=compartment)
    print("Prediction compelete")

    # Optionally squeeze the output
    if hyperparams.get('squeeze') is None or not hyperparams.get('squeeze'):
        img_arr = np.squeeze(img_arr)
        
    # save the output as a tiff
    # outfile = os.path.join(arg_dict['output_directory'], arg_dict['output_name'])
    # tifffile.imwrite(outfile, output)
    output = args.output_directory
    id_code = args.id_code
    expand_nuclear_area = eval(args.expand)
    if expand_nuclear_area is not None and expand_nuclear_area != 0:
        img_arr = skimage.segmentation.expand_labels(img_arr, distance=expand_nuclear_area)

    #Save as .tif file
    tifffile.imwrite(f'{output}/segments_mesmer-{id_code}.tif', img_arr)

    #Calculate and save areas
    (unique, counts) = np.unique(img_arr, return_counts=True)
    areas = np.asarray((unique, counts)).T
    np.savetxt(f'{output}/areas_mesmer-{id_code}.csv', areas, delimiter=",")