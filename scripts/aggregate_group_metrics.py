import pandas as pd
import numpy as np
import anndata as ad
import txsim as tx
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate metrics across different runs')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Output data directory- ')
    
    
    args = parser.parse_args()
    
    data = args.data

    if not os.path.exists( os.path.join(data, "aggregated") ):
        os.makedirs( os.path.join(data, "aggregated") )

    replicate_folders = [f.path for f in os.scandir(data) if f.is_dir() and 'replicate' in f.path]
    metric_list = []
    name_list = []
    
    #Read all replicates and aggreated 
    for replicate in replicate_folders:
        metric_file_name = os.path.join(replicate, f"group_metrics.csv")
        metric_list.append( pd.read_csv(metric_file_name, index_col = [0,1]) )
        metric_list[-1].columns = [replicate.split("/")[-1] +"-"+ x for x in metric_list[-1].columns] #rename columns to include replicate number
        name_list.extend(metric_list[-1].columns)
    aggregated_metric = pd.read_csv(os.path.join(data, f"aggregated/group_metrics.csv"), index_col=[0,1])
    
    metric = tx.metrics.aggregate_group_metrics(metric_list, aggregated_metric, name_list=name_list)
    metric.to_csv(os.path.join(data, f"aggregated/aggregated_group_metrics.csv"))
    