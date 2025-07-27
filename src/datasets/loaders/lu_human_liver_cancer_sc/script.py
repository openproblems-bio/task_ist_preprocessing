from pathlib import Path 
from scipy.sparse import csr_matrix, hstack
import pandas as pd
import anndata as ad
import urllib.request
from datetime import datetime

## VIASH START
par = {
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2022Lu_human_liver_cancer_sc.h5ad",
    "dataset_id": "2022Lu_human_liver_cancer_sc",
    "dataset_name": "2022Lu_human_liver_cancer_sc",
    "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614",
    "dataset_reference": "https://doi.org/10.1038/s41467-022-32283-3",
    "dataset_summary": "This dataset contains scRNA-seq data from human liver cancer cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human liver cancer cells.",
    "dataset_organism": "Homo sapiens"
}
meta = {
    "temp_dir": "./temp/datasets/2022Lu_human_liver_cancer_sc",
}
## VIASH END

# helper variables
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_URLS = {
    "counts": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2EscRNAseq%2ES71915%2Ecount%2Etxt%2Egz",
    "obs": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2Emetadata%2Eupdated%2Etxt%2Egz"
}
FILE_PATHS = {
    "counts": TMP_DIR / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
    "obs": TMP_DIR / "GSE149614_HCC.metadata.updated.txt.gz",
}


# Download the data
print("Downloading data", flush=True)
for key, url in FILE_URLS.items():
    urllib.request.urlretrieve(url, FILE_PATHS[key])

# Read data
print("Reading count matrix", flush=True)
chunk_size = 200

genes = []
X_sparse = []

t0 = datetime.now()

for i,chunk in enumerate(pd.read_csv(FILE_PATHS["counts"], sep="\t", chunksize=chunk_size)):

    if i % 10 == 0:
        print("\t", datetime.now() - t0, " " , i, "/128")

    genes += chunk.index.tolist()
    
    if i == 0:
        X_sparse = csr_matrix(chunk.values.T)
        obs = chunk.columns.tolist() # it's the same in each chunk since all cells are loaded
    else:
        X_sparse = hstack([X_sparse, csr_matrix(chunk.values.T)])

    del chunk

print("Reading obs", flush=True)
df_obs = pd.read_csv(FILE_PATHS["obs"], sep="\t", index_col=0)

assert (obs == df_obs.index.tolist())

# Create adata
print("Creating adata", flush=True)
adata = ad.AnnData(
    X=None,
    obs=df_obs,
    var=pd.DataFrame(index=genes),
    layers={"counts": X_sparse}
)

# Rename fields
rename_obs_keys = { 
    "cell_type": "celltype", 
    "batch": "sample",
    "donor_id": "patient",
    "cancer_stage": "stage", # TODO: this is currently not in the config yaml - how to handle this?
    # #TODO "cell_type_unified" (?), maybe call the unified one "cell_type" and the original one "cell_type_level1"
    # other keys: "batch", "assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id"
    # "diseases_ontology_term_id", "is_primary_data", "organism_ontology_term_id", "self_reported_ethnicity", 
    # "self_reported_ethnicity_ontology_term_id", "sex_ontology_term_id", "suspension_type", 
    # "suspension_type_ontology_term_id", "tissue_ontology_term_id", "tissue_general_ontology_term_id", "soma_joinid"
}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})

# Add additional information to obs
#TODO: Finish up the terms according to the ontology
#Ontology schema currently (13.03.2025) used in openproblems (CELLxGENE schema v4.0.0): 
#https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md
#(mentioned here: https://openproblems.bio/documentation/reference/openproblems/src-datasets#file-format:-raw-dataset)
store_info = { 
    "dataset_id": "2022Lu_human_liver_cancer_sc", #"GSE176078",
    "assay": "Chromium Single-Cell v2 3’ Chemistry Library", # from metadata from GEO GSE176078 #TODO: ontology
    "tissue": "liver",
    "disease": "liver cancer",
    # "disease_ontology_term_id": "PATO:0000461", #TODO: ontology
    "organism": "Homo sapiens",
    # "organism_ontology_term_id": "NCBITaxon:10090", #TODO: ontology
    "tissue_general": "liver",
    # "tissue_general_ontology_term_id": "UBERON:0000955", #TODO: ontology
    "development_stage": "adult", 
    # "development_stage_ontology_term_id": "MmusDv:0000110" #TODO: ontology
}
for key, value in store_info.items():
    adata.obs[key] = pd.Categorical([value] * adata.n_obs, categories=[value])

# Remove undesired columns
for key in adata.obs.columns:
    if (key not in rename_obs_keys.keys()) and (key not in store_info.keys()):
        print(f"Removing .obs['{key}']")
        del adata.obs[key]

# Var
adata.var["gene_symbol"] = adata.var_names
# TODO: can we also get ensembl ids? (adata.var["feature_id"]) 

# Uns
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

# Delete files if requested
if not par["keep_files"]:
    print("Removing files", flush=True)
    for file in FILE_PATHS.values():
        if file.exists():
            print("\t...", file, flush=True)
            file.unlink()

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
