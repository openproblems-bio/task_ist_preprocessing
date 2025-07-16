from pathlib import Path
import os
import pandas as pd
import anndata as ad


## VIASH START

par = {
    "input": "ftp://anonymous@ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad",
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2024Zuani_human_nsclc_sc.h5ad",
    "dataset_id": "2024Zuani_human_nsclc_sc", 
    "dataset_name": "2024Zuani_human_nsclc_sc", 
    "dataset_url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526",
    "dataset_reference": "https://doi.org/10.1038/s41467-024-48700-8", 
    "dataset_summary": "This dataset contains scRNA-seq data from human lung cancer cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human lung cancer cells.",
    "dataset_organism": "Homo sapiens"
}

meta = {
    "temp_dir": "./temp/datasets/2024Zuani_human_nsclc_sc",
}

## VIASH END


# Helper variables
TMP_DIR = Path(meta["temp_dir"] or "./tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
#FILE_PATH = TMP_DIR / par["input"].split("/")[-1]
#DOWNLOAD_URL = par["input"]
# Hard code the download url since it fails running on seqera: (TODO: either keep it like this and remove the input argument or remove this and resolve the issue)
#   (--2025-07-16 05:34:27--  http://_viash_par/input_1/10X_Lung_Tumour_Annotated_v2.h5ad
#   Resolving _viash_par (_viash_par)... failed: Name or service not known.
#   wget: unable to resolve host address ‘_viash_par’)
FILE_PATH = TMP_DIR / "10X_Lung_Tumour_Annotated_v2.h5ad"
DOWNLOAD_URL = "ftp://anonymous@ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad"

# Download the data (55GB)
os.system(f'wget "{DOWNLOAD_URL}" -P "{TMP_DIR}/"')
# os.system(f'wget "{DOWNLOAD_URL}" -P "{TMP_DIR}/" --show-progress')
adata = ad.read_h5ad(FILE_PATH)
# adata = adata[::100]

# Filter genes (not needed)
# sc.pp.filter_genes(adata, min_counts=1)

# Filter cells to NSCLC (~200k cells filtered out)
tumour_type_to_nsclc_status = {
    "NSCLC": "NSCLC",
    "Squamous cell carcinoma": "NSCLC",
    "Squamous dysplasia": "not NSCLC",  # pre-malignant lesion
    "Squamous cancer": "NSCLC", 
    "Adenocarcinoma ": "NSCLC",
    "Adenocarcinoma": "NSCLC",
    "NA": "not NSCLC",  # unclear / missing data
    "Mucinouse\nadenocarcinoma": "NSCLC",  
    "Presumed Lung cancer": "not NSCLC",  # not a confirmed subtype
    "Squamous carcinoma": "NSCLC",
    "Squamous cell lung cancer": "NSCLC",
    "lung adenocarcinoma": "NSCLC",
    "TTF1 +ve lung adenocarcinoma": "NSCLC",
    "Lung cancer": "not NSCLC"  # too generic #TODO: check from paper if this refers to NSCLC or not
}
adata.obs["NSCLC"] = adata.obs["tumour type"].map(tumour_type_to_nsclc_status)
adata = adata[adata.obs["NSCLC"] == "NSCLC"]

# Filter out cell types that should be removed
to_remove = adata.obs["Cell types"].str.endswith("(to remove)")
adata = adata[~to_remove]


# Rename or copy obs columns
rename_obs_keys = {
    "cell_type": 'Cell types',
    "donor_id": "patient",
    "sex": "sex",
    "batch": "batch",
}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})

# Store obs metadata with single values
store_info = { 
    "dataset_id": par["dataset_id"],
    "tissue": "lung",
    "disease": "NSCLC", 
    "organism": "Homo sapiens",
    "tissue_general": "lung",
    "development_stage": "adult", 
    # #TODO other keys: "assay", "assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id"
    # "diseases_ontology_term_id", "is_primary_data", "organism_ontology_term_id", "self_reported_ethnicity", 
    # "self_reported_ethnicity_ontology_term_id", "sex_ontology_term_id", "suspension_type", 
    # "suspension_type_ontology_term_id", "tissue_ontology_term_id", "tissue_general_ontology_term_id", "soma_joinid"
}
for key, value in store_info.items():
    adata.obs[key] = pd.Categorical([value] * adata.n_obs, categories=[value])

# Subset obs columns
obs_cols = list(rename_obs_keys.keys()) + list(store_info.keys())
adata.obs = adata.obs[obs_cols]

# Save uns metadata
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

# Add gene symbol column
adata.var["gene_symbol"] = adata.var_names

# Subset var columns
var_cols = ["gene_symbol"]
adata.var = adata.var[var_cols]

# Add layers
adata.layers['counts'] =  adata.X
del adata.X

# Delete files if requested
if not par["keep_files"]:
    print("Removing files", flush=True)
    if FILE_PATH.exists():
        print("\t...", FILE_PATH, flush=True)
        FILE_PATH.unlink()

# Write adata
print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
