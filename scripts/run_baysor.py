import pandas as pd
import argparse
import os
from pathlib import Path


def convert_str_ids_to_ints(df, file_path_for_error_messages=None):
    """Convert cell ids like "CR4b68f93d8-27" to 27
    
    The file argument is just for creating more informative Error messages.
    """
    
    df = df.copy()
    file = file_path_for_error_messages
    
    unique_cell_values = df.loc[~df["cell"].isnull(), "cell"].unique()
    unique_types = {type(value) for value in unique_cell_values}
    n_cells_pre = len(df.loc[~df["cell"].isnull(),"cell"].unique())
    if (len(unique_types) == 1) and (str in unique_types):
        df.loc[~df["cell"].isnull(),"cell"] = df.loc[~df["cell"].isnull(),"cell"].apply(lambda i: i.split("-")[-1]).astype(int)
    elif (len(unique_types) != 1):
        raise ValueError(f"Non NaN values of column 'cell' in file {file} have multiple types: {unique_types}")
    n_cells_post = len(df.loc[~df["cell"].isnull(),"cell"].unique())
    
    if n_cells_pre != n_cells_post:
        raise ValueError(f"Number of cells changed after conversion to integers, probably baysor used different substrings with same integers for some cells. Check file: {file}")
        
    # Convert nan values to 0 (background)
    df.loc[df["cell"].isnull(),"cell"] = 0
    
    # Conver nan vlaues to None for the cell type column #TODO: check why we need this here for baysor
    if "celltype" in df.columns:
        df.loc[df["celltype"].isnull(),"celltype"] = "None_sp"
        
    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells using Baysor')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]') # spot.csv
    parser.add_argument('-o', '--output', required=False, tpye=str, default=None, help="Assignment output csv path")
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image') # <results>/<dataset>
    parser.add_argument('-s', '--segment', default=None, type=str,
        help='Segmentation method used for image') # DOES THIS CONTAIN ID? I think yes
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') # OK
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving') # clear
    parser.add_argument('--temp', default=None, type=str, 
        help='Temp data directory for intermediate files') # 
    
    DEFAULT_HYPERPARAMS = {
        # Taken from https://github.com/kharchenkolab/Baysor/blob/master/configs/example_config.toml
        
        # [Data]
        # Name of the x column in the input data. Default: "x"
        #'x-column' : '"x"',
        'x' : '"x"',
        # Name of the y column in the input data. Default: "y"
        #"y-column" : '"y"',
        "y" : '"y"',
        # Name of the y column in the input data. Default: "z"
        #"z-column" : '"z"',
        "z" : '"z"',
        # Name of gene column in the input data. Default: "gene"
        #"gene-column" : '"Gene"',
        "gene" : '"Gene"',
        # Minimal number of molecules per gene. Default: 1
        "min_molecules_per_gene" : 1,
        # Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Default: 3
        "min_molecules_per_cell" : 3,
        # Scale parameter, which suggest approximate cell radius for the algorithm. This parameter is required.
        "scale" : 50,
        # Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with "%" to set it relative to scale. Default: "25%"
        # "scale-std" : '"25%"',
        
        # Not exactly sure if this one should be in [Data], therefore we don't provide it via the toml, so not possibl issues here.
        "prior-segmentation-confidence" : 0.2,
        
        # # Use scale estimate from DAPI if provided. Default: true
        # "estimate-scale-from-centers" : "true",
        # # Minimal number of molecules in a segmented region, required for this region to be considered as a possible cell. Default: min-molecules-per-cell / 4
        # # "min-molecules-per-segment" : 2, # TODO: throws Error ErrorException("Unexpected value in the config: 'min-molecules-per-segment'")
        # 
        # #[Sampling] These parameters shouldn't normally be changed
        # # Prior weight of assignment a molecule to new component. Default: 0.2
        # "new-component-weight" : 0.2,
        # # Fraction of distributions, sampled at each stage. Default: 0.3
        # "new-component-fraction" : 0.3,
        
        #[Plotting]
        # Number of neighbors (i.e. 'k' in k-NN), which is used for gene composition visualization. Larger numbers leads to more global patterns. Default: estimate from min-molecules-per-cell
        # "gene-composition-neigborhood" : 20, # TODO: ErrorException("Unexpected value in the config: 'gene-composition-neigborhood'")
        # Number of pixels per cell of minimal size, used to estimate size of the final plot. For most protocols values around 7-30 give enough visualization quality. Default: 15
        # "min-pixels-per-cell" : 15, # TODO: didn't check this one, plotting is not so important i guess
    }

    
    
    
    
    args = parser.parse_args()

    molecules = args.molecules
    assign_out = args.output
    data = args.data
    segmentation_method = args.segment
    hyperparams = eval(args.hyperparams)
    hparams = {}
    for key in DEFAULT_HYPERPARAMS:
        if hyperparams is not None and key in hyperparams:
            hparams[key] = hyperparams[key]
        else:
            hparams[key] = DEFAULT_HYPERPARAMS[key]
    id_code = args.id_code
    segment = True if args.segment is not None else False
    temp = args.temp if args.temp is not None else data
    temp = Path(temp)
    toml_file = temp / 'config.toml'
    #temp = Path(temp) / f"baysor_{id_code}"
    #toml_file = temp / 'config.toml'

    temp.mkdir(parents=True, exist_ok=True)
    #if not os.path.exists(temp):
    #    os.makedirs(temp)

    baysor_seg = os.path.join(temp, "segmentation.csv")
    baysor_cell = os.path.join(temp, "segmentation_cell_stats.csv")
    
    # Remove existing outputs (otherwise Errors while running Baysor might be overseen)
    if os.path.isfile(baysor_seg):
        os.remove(baysor_seg)
    if os.path.isfile(baysor_cell):
        os.remove(baysor_cell)
        
    # Write toml
    with open(toml_file, "w") as file:
        for key, val in hparams.items():
            # TODO: extend the toml headers, check https://github.com/LouisK92/Baysor/blob/master/configs/example_config.toml
            #if key == "x-column":
            if key == "x":
                #file.write(f'[Data]\n')
                file.write(f'[data]\n')
            elif key == "new_component_weight":
                #file.write(f'\n[Sampling]\n')
                file.write(f'\n[sampling]\n')
            if key not in ["scale", "prior-segmentation-confidence"]:
                file.write(f'{key} = {val}\n')
    
    # Note: we provide scale separately because when providing it via .toml baysor can complain that's it's not a float
    baysor_cli = f"run -s {hparams['scale']} -c {toml_file} -o {temp}/ {molecules}"
    #baysor_cli = "run "
    #for key, val in hyperparams.items():
    #    baysor_cli += f"--{key} {val} "
        
    if segment:
        print("Running Baysor with prior segmentation")
        baysor_cli += f" --prior-segmentation-confidence {hparams['prior-segmentation-confidence']}"
        baysor_cli += f" {data}/segments_{segmentation_method}.tif"
        #baysor_cli += f" {data}/segments_{segmentation_method}.ome.tif" #TODO: ideally the container works with this line
        #baysor_cli += f"-o {temp}/ {molecules} {data}/segments_{segmentation_method}.tif"
        #baysor_cli += f"{molecules} -o {temp} --save-polygons=geojson -p {data}/segments_{segmentation_method}.tif"
            
    else:
        print("Running Baysor without prior segmentation")
        #baysor_cli += f"-o {temp}/ {molecules}"
        #baysor_cli += f"{molecules} -o {temp} --save-polygons=geojson"

    #os.system(f'''/Baysor/bin/baysor {baysor_cli}''') 
    print(f'''baysor {baysor_cli}''')
    os.system(f'''baysor {baysor_cli}''')
    print("Ran Baysor")

    df = pd.read_csv(baysor_seg)
    spots = pd.read_csv(molecules)
    spots["cell"] = df["cell"]
    
    spots.rename(columns = {spots.columns[0]:'Gene', spots.columns[1]:'x', spots.columns[2]:'y'} , inplace = True)

    df = pd.read_csv(baysor_cell)
    areas = df[['cell','area']]

    spots = convert_str_ids_to_ints(spots, file_path_for_error_messages=baysor_seg)
    areas = convert_str_ids_to_ints(areas, file_path_for_error_messages=baysor_cell)    
    
    #Save to csv (TODO: added the -o input argument, should probably only use that one, not allowing "None")
    if assign_out is None:
        if segment:
            areas.to_csv(f'{data}/areas_{segmentation_method}_baysor-{id_code}.csv', index = False, header = False)
            spots.to_csv(f'{data}/assignments_{segmentation_method}_baysor-{id_code}.csv', index = False)
        else:
            areas.to_csv(f'{data}/areas_baysor-{id_code}.csv', index = False, header = False)
            spots.to_csv(f'{data}/assignments_baysor-{id_code}.csv', index = False)
    else:
        spots.to_csv(assign_out, index = False)
        areas_out = assign_out.replace("assignments_","areas_")
        areas.to_csv(areas_out, index = False, header = False)
        
