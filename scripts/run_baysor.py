import pandas as pd
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells using Baysor')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-s', '--segment', default=None, type=str,
        help='Segmentation method used for image')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    parser.add_argument('--temp', default=None, type=str, 
        help='Temp data directory for intermediate files')
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    segmentation_method = args.segment
    hyperparams = eval(args.hyperparams)
    id_code = args.id_code
    segment = True if args.segment is not None else False
    temp = args.temp if args.temp is not None else data

    if not os.path.exists(temp):
        os.makedirs(temp)

    os.system('''julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="PackageCompiler"))'
    julia -e 'using Pkg; Pkg.add([PackageSpec(name="UMAP", rev="master"), PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")])'
    julia -e 'import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); Pkg.instantiate();'
    julia -e 'using PackageCompiler; import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor))))'
    ''')

    if segment:
        print("Running Baysor with prior segmentation")
        baysor_cli = f"run -s 40 {molecules} -o {temp} " +
            f"--save-polygons=geojson -p {data}/segments_{segmentation_method}.tif"
            
    else:
        print("Running Baysor without prior segmentation")
        baysor_cli = f"run -s 40 {molecules} -o {temp} " +
            "--save-polygons=geojson -p "
    os.system(f'''julia -E 'import Baysor; Baysor.run_cli("{baysor_cli}")' ''')
    print("Ran Baysor")

    baysor_seg = os.path.join(temp, "segmentation.csv")
    baysor_cell = os.path.join(temp, "segmentation_cell_stats.csv")

    df = pd.read_csv(baysor_seg)
    assignments = df[['gene','x','y','cell']]
    assignments.columns = ["Gene","x","y","cell"]

    df = pd.read_csv(baysor_cell)
    areas = df[['cell','area']]

    #Save to csv
    if segment:
        areas.to_csv(f'{data}/areas_{segmentation_method}_baysor-{id_code}.csv', index = False, header = False)
        assignments.to_csv(f'{data}/assignments_{segmentation_method}_baysor-{id_code}.csv', index = False)
    else:
        areas.to_csv(f'{data}/areas_baysor-{id_code}.csv', index = False, header = False)
        assignments.to_csv(f'{data}/assignments_baysor-{id_code}.csv', index = False)

    #cell_types.to_csv(f'{data}/celltypes_{segmentation_method}_pciSeq-{id_code}.csv', header = False)
    
    