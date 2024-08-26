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
    parser.add_argument('-m', '--methods', required=True, type=str,
        help='Methods used to generate count matrices')
    parser.add_argument('-t', '--metric_type', required=True, type=str,
        help='Type of metrics, either "quality_metrics" or "metrics" ')
    
    
    args = parser.parse_args()
    
    data = args.data
    methods = args.methods
    metric_type = args.metric_type

    if not os.path.exists( os.path.join(data, "aggregated") ):
        os.makedirs( os.path.join(data, "aggregated") )

    replicate_folders = [f.path for f in os.scandir(data) if f.is_dir() and 'replicate' in f.path]
    metric_list = []
    name_list = []
    
    for replicate in replicate_folders:
        metric_file_name = os.path.join(replicate, f"{metric_type}_{methods}.csv")
        metric_list.append( pd.read_csv(metric_file_name, index_col = 0) )
        name_list.append(replicate.split("/")[-1])

    aggregated_metric = pd.read_csv(os.path.join(data, f"aggregated/{metric_type}_{methods}.csv"), index_col=0)
    metric = tx.metrics.aggregate_metrics(metric_list, aggregated_metric, name_list=name_list)
    metric.to_csv(os.path.join(data, f"aggregated/aggregated_{metric_type}_{methods}.csv"))
    
    #unfiltered metrics

    if not os.path.exists( os.path.join(data, "aggregated/unfiltered") ):
        os.makedirs( os.path.join(data, "aggregated/unfiltered") )
    
    metric_list = []
    name_list = []
    for replicate in replicate_folders:
        metric_file_name = os.path.join(replicate, f"unfiltered/{metric_type}_{methods}.csv")
        metric_list.append( pd.read_csv(metric_file_name, index_col = 0) )
        name_list.append(replicate.split("/")[-1])

    aggregated_metric = pd.read_csv(os.path.join(data, f"aggregated/unfiltered/{metric_type}_{methods}.csv"), index_col=0)
    metric = tx.metrics.aggregate_metrics(metric_list, aggregated_metric, name_list=name_list)
    metric.to_csv(os.path.join(data, f"aggregated/unfiltered/aggregated_{metric_type}_{methods}.csv"))
