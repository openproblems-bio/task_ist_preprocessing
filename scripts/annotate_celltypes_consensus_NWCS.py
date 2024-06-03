
from pathlib import Path
import argparse
import pandas as pd



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate cells with celltype labels using NWCS consensus')

    parser.add_argument('-i', '--input_files', nargs='+', help='List of CSV file paths to the outputs of the different methods')
    parser.add_argument('-o', '--output_file', default='output.csv', help='Output file path')
    #parser.add_argument('-p', '--hyperparams', default=None, type=str,help='Optional dictionary (as string) of method-specific parameters') 
    args = parser.parse_args()

    
    # Create output directory if it does not exist
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)

    # Load default hyperparameter
    # hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    # with open(hparams_defaults_csv, 'r',encoding='utf-8') as file:
    #     defaults = yaml.safe_load(file)
    #     hparams_defaults = defaults["nwconsensus"] 
    
    # # Parse parameters
    # hyperparams = eval(args.hyperparams) if args.hyperparams else {}
    # hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    # hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}



    df_dict = {file: pd.read_csv(file) for file in args.input_files}



    all_celltypes = set()
    for df in df_dict.values():
        all_celltypes.update(df['celltype'].unique())

    result_dfs = {}
    desired_columns = ['cell_id', 'celltype', 'score']
    for method, df in df_dict.items():
        df = df[desired_columns]
        # Sort cell_id to ensure consistent order across all resulting DataFrames
        df_sorted = df.sort_values('cell_id')
        # Create a pivot table with cell_id as row index, celltype as columns, and score as values
        pivot_df = df_sorted.pivot(index='cell_id', columns='celltype', values='score')        
         # Reindex to ensure all cell types appear in the same order across all resulting DataFrames
        pivot_df = pivot_df.reindex(columns=all_celltypes)
        # Fill missing values with -1
        pivot_df.fillna(-1, inplace=True)
        result_dfs[method] = pivot_df

  
    #calculate average
    matrix_sum = sum(df.values for df in result_dfs.values())
    average_matrix = matrix_sum / len(result_dfs)
    average_df = pd.DataFrame(average_matrix, index=result_dfs[next(iter(result_dfs))].index, columns=result_dfs[next(iter(result_dfs))].columns)


    # Find the celltype with the maximum score for each cell_id
    max_celltype = average_df.idxmax(axis=1)
    consensus_df = max_celltype.to_frame().reset_index()
    consensus_df.columns = ['cell_id', 'celltype']
    consensus_df['score'] = average_df.max(axis=1).values

    

      
    # Save annotation
    consensus_df.to_csv(args.output_file, index=False)
