import argparse
import tifffile
import pandas as pd
import txsim as tx



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Subset data to tile')
    parser.add_argument('-i', '--image', required=True, type=str, help='Input ome tif') 
    parser.add_argument('-m', '--molecules', required=False, default=None, type=str, help='Input spots csv') 
    parser.add_argument('-t', '--tile_info', required=True, type=str,
        help='Numbers that define the tile given as {ny}-{nx}-{y_id}-{x_id}-{n_expand_px}') 
    parser.add_argument('-g', '--output_img', required=True, type=str, help='Output image')
    parser.add_argument('-o', '--output_mol', required=False, default=None, type=str, help='Output molecules')
    
    # Get arguments
    args = parser.parse_args()
    
    img_path = args.image
    spots_path = args.molecules
    ny, nx, y, x, n_expand_px = [int(s) for s in args.tile_info.split("-")]
    img_out_path = args.output_img
    spots_out_path = args.output_mol
    
    assert (spots_path is None) or ((spots_path is not None) and (spots_out_path is not None)); "If a spots table is given a spots output path also needs to be given."
    
    # Get tile intervals
    y_len, x_len = tx.data_utils.get_image_shape(img_path)
    x_min, x_max, y_min, y_max, offset_x, offset_y = tx.data_utils.get_tile_intervals(
        x, y, nx, ny, x_len, y_len, n_expand_px
    )
    
    # Load and subset data
    img = tifffile.imread(img_path)[y_min:y_max,x_min:x_max]
    if spots_path:
        spots = pd.read_csv(spots_path)
        spots = tx.data_utils.get_spots_tile(spots, x_min, x_max, y_min, y_max)
    
    # Save data
    with tifffile.TiffWriter(img_out_path, bigtiff=True) as tif:
        # TODO: as we have an ome tif as input we can just use the same xml meta data! 
        #       currently not relevant, but will be when converting to physical sizes.
        metadata={ 
            'PhysicalSizeX': 1,
            'PhysicalSizeXUnit': 'um',
            'PhysicalSizeY': 1,
            'PhysicalSizeYUnit': 'um'
        }
        tif.write(
            img,
            metadata=metadata
        )
    tifffile.imwrite(img_out_path.replace("ome.tif","tif"), img)
    if spots_path:
        spots.to_csv(spots_out_path, ignore_index=True)
    
