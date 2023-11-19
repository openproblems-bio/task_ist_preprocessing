import argparse
import numpy as np
import pandas as pd
import txsim as tx


def extract_tile_info_from_path(path):
    """extract info from strs like '.../assignments_{method}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    """
    strs = path.rsplit("_",5)
    ny = int(strs[1][2:])
    nx = int(strs[2][2:])
    y_id = int(str[3])
    x_id = int(str[4])
    n_expand_px = int(str[5].split(".")[0])
    
    return ny, nx, y_id, x_id, n_expand_px


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Aggregate assignment outputs of different tiles')
    parser.add_argument('-p', '--input', nargs='+', required=True, 
                        help="Assignments (and areas) input file paths of tiles")
    parser.add_argument('-o', '--output', nargs='+', required=True, 
                        help="Assignments (and areas) output file paths for aggregated csv")
    parser.add_argument('-m', '--molecules', required=True, help="Spots paths (needed to keep the spots order)")
    parser.add_argument('-i', '--image', required=True, help="Image path (to get img dimensions)")
    
    # Get arguments
    args = parser.parse_args()
    
    assign_paths = [f for f in args.input if ("assignments_" in f)]
    area_paths = [f for f in args.input if ("areas_" in f)]
    spots_path = args.molecules
    img_path = args.image
    out_paths = args.output
    if type(out_paths) == str:
        assert "assignments" in out_paths; f"If only 1 path is given ({out_paths}), it should be the assignments path"
        spots_out_path = out_paths
        areas_out_path = None
    elif (type(out_paths) == list) and (len(out_paths) == 2):
        spots_out_path = [p for p in out_paths if "assignments" in p][0]
        areas_out_path = [p for p in out_paths if "areas" in p][0]
        assert len(assign_paths) == len(area_paths)
    else:
        raise ValueError(f"Unexpected output paths: {out_paths}")
    
    # Load spots and image shape
    spots = pd.read_csv(spots_path)
    y_len, x_len = tx.data_utils.get_image_shape(img_path)
    
    # Add spots columns (tile_id will be removed later)
    spots["cell"] = None
    spots["not_assigned_due_to_tiling"] = False
    spots["tile_id"] = None
    
    area_path = None
    area_df = None
    areas = pd.DataFrame(columns=["area","tile_id"])
    
    for assign_path in assign_paths:
        ny, nx, y, x, n_expand_px = extract_tile_info_from_path(assign_path)
        if areas_out_path is not None:
            area_path = [p for p in area_paths if p.endswith(f"_ny{ny}_nx{nx}_{y}_{x}_px{n_expand_px}.csv")][0]
        
        # Load
        tile_df = pd.read_csv(assign_path)
        if area_path:
            area_df = pd.read_csv(area_path, index_col=0, header=None)
            area_df.index = area_df.index.astype(int)
            area_df.columns = ["area"]
        
        # Translate back the spots coordinates
        x_min_w_exp, _, y_min_w_exp, _, _, _ = tx.data_utils.get_tile_intervals(x, y, nx, ny, x_len, y_len, n_expand_px)
        tile_df["x"] += x_min_w_exp
        tile_df["y"] += y_min_w_exp
        
        # get list of cell ids that cross the borders, set those as background, and renumber cells
        x_min, x_max, y_min, y_max, _, _ = tx.data_utils.get_tile_intervals(x, y, nx, ny, x_len, y_len, 0)
        tile_df["in_tile"] = tx.data_utils.get_tile_mask(tile_df, x_min, x_max, y_min, y_max)
        crossing_cells = tile_df.loc[~tile_df["in_tile"],"cell"].unique().tolist()
        if 0 in crossing_cells:
            crossing_cells.remove(0)
        tile_df["not_assigned_due_to_tiling"] = tile_df['cell'].isin(crossing_cells)
        # set unassigned as background
        tile_df.loc[tile_df["not_assigned_due_to_tiling"],'cell'] = 0
        # renumber cell ids
        ids = sorted(tile_df["cell"].unique())
        tile_df["cell"] = tile_df["cell"].apply(lambda x: ids.index(x))
        if isinstance(area_df, pd.DataFrame):
            area_df = area_df.loc[~area_df.index.isin(crossing_cells)]
            area_df.index = area_df.index.to_series().apply(lambda x: ids.index(x))
        
        # Subset to spots in tile
        tile_df = tile_df.loc["in_tile"]
        del tile_df["in_tile"]
        
        # Add cell ids (and other columns) to spots
        mask = tx.data_utils.get_tile_mask(spots, x_min, x_max, y_min, y_max)
        assert mask.sum() == len(tile_df)
        assert np.all(spots.loc[mask,"Gene"] == tile_df["Gene"])
        spots.loc[mask,"cell"] = tile_df["cell"].values
        spots.loc[mask,"not_assigned_due_to_tiling"] = tile_df["not_assigned_due_to_tiling"].values
        # Set a tile_id since we'll need to renumber cells per tile when tiles are aggregated
        spots.loc[mask,"tile_id"] = y*nx + x
        if isinstance(area_df, pd.DataFrame):
            area_df["tile_id"] = y*nx + x
            areas = pd.concat([areas,area_df])
        
    # Renumber cell ids
    max_id = 0
    for tile_id in range(nx*ny):
        spots.loc[(spots["tile_id"]==tile_id) & (spots["cell"]!=0),"cell"] += max_id
        if areas_out_path is not None:
            areas.loc[(areas["tile_id"]==tile_id) & (areas.index!=0),"cell"] += max_id
        max_id = spots.loc[(spots["tile_id"]==tile_id) & (spots["cell"]!=0),"cell"].max()
        
    # Delete tile_id
    del spots["tile_id"]
        
    # Save assignments
    spots.to_csv(spots_out_path, ignore_index=True)
    # Save areas
    if areas_out_path is not None:
        if (0 in areas.index) and (np.sum(areas.index == 0) > 1):
            # Check if background (0) occurs multiple times. If yes, sum up #TODO: I think backgrounds should be filtered out in previous steps
            areas = pd.concat([
                pd.DataFrame(index=[0],data={"area":[areas.loc[0,"area"].sum()]}), 
                areas.loc[areas.index != 0]
            ])
        areas["area"].to_csv(areas_out_path)
        
