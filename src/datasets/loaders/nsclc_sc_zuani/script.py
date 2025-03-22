from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
import anndata as ad
import os
import sys
import scipy.sparse
import scanpy

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
    "temp_dir": "./tmp/nsclc_sc_zuani/",
}

## VIASH END


#os.system(f'wget "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad" -P {meta["temp_dir"]}')
#adata = ad.read_h5ad( f'{meta["temp_dir"]}10X_Lung_Tumour_Annotated_v2.h5ad', backed="r")


TMP_DIR = Path(meta["temp_dir"] or "./tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_PATHS = {"file": TMP_DIR / "cropped_sc.h5ad"}
os.system(f'wget http://192.168.2.46:8000/file/cropped_sc.h5ad -P ./tmp/')
adata = ad.read_h5ad( './tmp/cropped_sc.h5ad')

genes_sum = adata.X.toarray().sum(0)
adata = adata[:, genes_sum != 0]

rename_obs_keys = {
    "cell_type": 'Cell types'}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})



store_info = { 
    "dataset_id": "E-MTAB-13526",
    "assay": "Unknown", 
    "sex": "female, male", 
    "tissue": "lung",
    "disease": "lung adenocarcinoma, normal, non-small cell lung cancer, lung squamous cell carcinoma",
    "organism": "Homo sapiens",
    "tissue_general": "lung",
    "development_stage": "adult", 
}
for key in ["dataset_id", "tissue", "organism", "tissue_general", "development_stage"]:
    adata.obs[key] = pd.Categorical([store_info[key]] * adata.n_obs, categories=[store_info[key]])

uns_info = { "dataset_id": "E-MTAB-13526" ,
              "dataset_name":"E-MTAB-13526" , 
              "dataset_url":"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526" ,
              "dataset_reference": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526",
              "dataset_summary": 'none',
              "dataset_description":'none',
              "dataset_organism": 'Homo sapiens' 
}

for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = uns_info[key]

adata.var["gene_symbol"] = adata.var_names

adata.layers['counts'] =  adata.X.toarray()
del adata.X

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
