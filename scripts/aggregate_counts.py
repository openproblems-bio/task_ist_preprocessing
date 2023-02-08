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
    parser.add_argument('-f', '--files', required=True, type=str,
        help='List of all files or "all" for all runs in dataset')
    
    args = parser.parse_args()
    
    data = args.data
    files = args.files


    