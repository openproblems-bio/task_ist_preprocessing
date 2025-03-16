from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
import anndata as ad
import os
import sys

## VIASH START

par = {
    "output": "EMTAB13526_Lung_sc.h5ad",
    "dataset_id": "E-MTAB-13526",
    "dataset_name": "E-MTAB-13526",
    "dataset_url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526",
    "dataset_reference": "https://doi.org/10.1038/s41467-024-48700-8", 
    "dataset_summary": "This dataset contains scRNA-seq data from human lung cancer cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human lung cancer cells.",
    "dataset_organism": "Homo sapiens"
}

meta = {
    "temp_dir": "./tmp/nsclc_sc_zuani",
}

## VIASH END


#os.system(f'wget "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad" -P {par["output"]}')


adata = ad.read_h5ad( '/mnt/d/file_h5ad/file.h5ad/10X_Lung_Tumour_Annotated_v2.h5ad', backed="r")
adata = adata[:100].to_memory()

rename_obs_keys = {
    "cell_type": 'Cell types'}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})



print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
