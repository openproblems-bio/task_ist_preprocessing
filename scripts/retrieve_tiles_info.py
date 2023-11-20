
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import txsim as tx


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Aggregate assignment outputs of different tiles')
    parser.add_argument('-o', '--output', required=True, 
                        help="Output path for tile info csv.")
    parser.add_argument('-m', '--molecules', required=True, help="Spots paths (needed to keep the spots order)")
    parser.add_argument('-i', '--image', required=True, help="Image path (to get img dimensions)")
    parser.add_argument('-nb', '--n_spots_per_tile_baysor', required=True, type=int, 
                        help="Maximal number of spots per tile (baysor)")
    parser.add_argument('-nc', '--n_spots_per_tile_clustermap', required=True, type=int, 
                        help="Maximal number of spots per tile (clustermap)")
    parser.add_argument('-ne', '--n_pixel_extend_tile', required=True, type=int, 
                        help="Number of pixels added at tile borders to account for edge effects")
    
    
    # Get arguments
    args = parser.parse_args()
    
    spots_path = args.molecules
    img_path = args.image
    out_path = args.output
    n_baysor = args.n_spots_per_tile_baysor
    n_clustermap = args.n_spots_per_tile_clustermap
    extend_n_px = args.n_pixel_extend_tile
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Init info df
    df = pd.DataFrame(index=["baysor","clustermap"], columns=["ny","nx","extend_n_px"])
    df["extend_n_px"] = extend_n_px
    
    # Load input data
    img_shape = tx.data_utils.get_image_shape(img_path)
    spots = pd.read_csv(spots_path)

    # Find optimal tiles
    ny, nx = tx.data_utils.find_optimal_tile_division_for_nspots_limit(
        img_shape, spots, n_spots_max=n_baysor, extend_n_pixels=extend_n_px
    )
    df.loc["baysor","ny"] = ny
    df.loc["baysor","nx"] = nx
    
    if n_clustermap != n_baysor:
        ny, nx = tx.data_utils.find_optimal_tile_division_for_nspots_limit(
            img_shape, spots, n_spots_max=n_baysor, extend_n_pixels=extend_n_px
        )
    df.loc["clustermap","ny"] = ny
    df.loc["clustermap","nx"] = nx

    # Save
    df.to_csv(out_path)
    
