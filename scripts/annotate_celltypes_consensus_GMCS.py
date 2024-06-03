from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import os



def geometric_median(points, max_iter=100, tol=1e-5):
    """
    Compute the geometric median of a set of points using Weiszfeld's algorithm.
    
    Args:
        points (numpy.ndarray): Array of shape (m, n) representing m points in n-dimensional space.
        max_iter (int): Maximum number of iterations for the algorithm.
        tol (float): Tolerance for stopping criterion.
        
    Returns:
        numpy.ndarray: Geometric median point in n-dimensional space.
    """
    num_points, dim = points.shape
    weights = np.ones(num_points) / num_points
    prev_median = np.zeros(dim)
    
    for _ in range(max_iter):
      
        median = np.sum(weights[:, None] * points, axis=0) / np.sum(weights)
        
       
        distances = np.linalg.norm(points - median, axis=1)
        weights = 1 / (distances + np.finfo(float).eps)
        
        if np.linalg.norm(prev_median - median) < tol:
            break
        
        prev_median = median
    
    return median



def read_data(file_path):
    """
    Read data from a CSV file and return the DataFrame without 'celltype' and 'score' columns, add columns for cell type probabilities if they don't exist
    
    Args:
        file_path (str): Path to the CSV file.
        
    Returns:
        pandas.DataFrame: DataFrame containing the data 
    """
    data = pd.read_csv(file_path)
    # Check if there are any columns apart from 'celltype', 'score', and 'cell_id'
    other_columns = [col for col in data.columns if col not in ['celltype', 'score', 'cell_id']]

    if not other_columns:
        # Create new columns for each unique value in 'celltype'
        unique_cell_types = data['celltype'].unique()
        for cell_type in unique_cell_types:
            # For rows where 'celltype' is equal to the current cell_type,
            # set the value of the corresponding column to 'score'
            data[cell_type] = data.loc[data['celltype'] == cell_type, 'score']

    # Drop 'celltype' and 'score' columns
    data.drop(columns=['celltype', 'score'], inplace=True, errors='ignore')

    return data

def geometric_median_combining(method_files):
    """
    Combine multiple methods using the geometric median combining strategy.
    
    Args:
        method_files (list of str): List of file paths to the CSV files containing method data.
        
    Returns:
        pandas.DataFrame: DataFrame containing cell_id and combined_matching.
    """




    # Initialize list to store cell_ids for each method
    cell_ids_list = []

    # Read data from each method file
    matchings = []
    all_cell_types = set()

    for file_path in method_files:
        print(file_path)
        method_data = read_data(file_path)
        # Sort method data by cell_id to ensure correct combination
        method_data.sort_values(by='cell_id', inplace=True)
        cell_ids_list.append(method_data['cell_id'])
        matchings.append(method_data)
        all_cell_types.update(method_data.columns)
        
        

    
    cell_ids = cell_ids_list[0] #we assume all files have the same cell_ids 
    all_cell_types_list = ['cell_id'] + list(all_cell_types - {'cell_id'})#make cell_id the first column
 
   
    for i, matching in enumerate(matchings):
        missing_cell_types = all_cell_types - set(matching.columns)

        # add missing cell types
        missing_df = pd.DataFrame({cell_type: [np.nan] * len(matching) for cell_type in missing_cell_types})
        matching = pd.concat([matching, missing_df], axis=1)

        # Reorder columns to match all_cell_types_list order
        matching = matching.reindex(columns=all_cell_types_list)


        # Fill nan values by  distributed remaining probability
        row_sums = matching.drop(columns='cell_id').sum(axis=1)
        remaining_proportion = 1 - row_sums
        n_missing_cell_type_columns = matching.isna().sum(axis=1)
        adjustment_values = remaining_proportion / n_missing_cell_type_columns
        adjustment_values[n_missing_cell_type_columns == 0] = 0
        adjustment_values[remaining_proportion <= 0] = 0
  
        columns_to_fill = matching.columns[1:]
        for column in columns_to_fill:
            matching[column] = matching[column].fillna(adjustment_values)

        matchings[i] = matching


    # Create a 3D array
    n_cells = len(cell_ids)
    n_methods = len(matchings)
    n_celltypes = len(all_cell_types_list)-1
    data_3d = np.empty((n_cells, n_methods, n_celltypes))
    for i, matching in enumerate(matchings):
        for j, cell_id in enumerate(cell_ids):
            data_3d[j, i, :] = matching.loc[matching['cell_id'] == cell_id].iloc[0].drop('cell_id').values



    
   
    # Compute geometric medians for each cell
    geometric_medians = []
    
    for i in range (data_3d.shape[0]):
        geometric_medians.append(geometric_median(data_3d[i]))

   

    combined_df = pd.DataFrame(np.array(geometric_medians), columns=all_cell_types_list[1:])
   
    # Add cell_id column
    combined_df.insert(0, 'cell_id', cell_ids.values)
    return combined_df









if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate cells with celltype labels using consensus')

    parser.add_argument('-i', '--input_files', nargs='+', help='List of CSV file paths to the outputs of the different methods')
    parser.add_argument('-o', '--output_file', default='output.csv', help='Output file path')
    #parser.add_argument('-p', '--hyperparams', default=None, type=str,help='Optional dictionary (as string) of method-specific parameters') 
    args = parser.parse_args()

    
    # Create output directory if it does not exist
    

    output_directory = Path(args.output_file).parent
    if not output_directory.exists():
        output_directory.mkdir(parents=True)

   
    # # Load default hyperparameter
    # hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    # with open(hparams_defaults_csv, 'r',encoding='utf-8') as file:
    #     defaults = yaml.safe_load(file)
    #     hparams_defaults = defaults["consensus"] 
    
    # # Parse parameters
    # hyperparams = eval(args.hyperparams) if args.hyperparams else {}
    # hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    # hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}



    combined_matching = geometric_median_combining(args.input_files)
    combined_matching.insert(1, 'celltype', combined_matching.iloc[:, 1:].idxmax(axis=1))
    combined_matching.insert(2, 'score', combined_matching.iloc[:, 2:].max(axis=1))
  
   
    # Save annotation
    combined_matching.to_csv(args.output_file, index=False)



    







    