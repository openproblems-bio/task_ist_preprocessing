from pathlib import Path 
import pandas as pd
import anndata as ad
import scanpy as sc
import urllib.request
import tarfile


## VIASH START
par = {
    "cancer_subtypes": ['HER2+', 'TNBC', 'ER+'],
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2021Wu_human_breast_cancer_sc.h5ad",
    "dataset_id": "2021Wu_human_breast_cancer_sc",
    "dataset_name": "2021Wu_human_breast_cancer_sc",
    "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
    "dataset_reference": "https://doi.org/10.1038/s41588-021-00911-1", #TODO: bibtex not doi, also adjust config.vsh.yaml
    "dataset_summary": "This dataset contains scRNA-seq data from human breast cancer cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human breast cancer cells.",
    "dataset_organism": "Homo sapiens"
}
meta = {
    "temp_dir": "./temp/datasets/2021Wu_human_breast_cancer_sc",
}
## VIASH END

# helper variables
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_PATHS = {
    "tar_gz": TMP_DIR / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
    "counts_mtx": TMP_DIR / "Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx",
    "barcodes": TMP_DIR / "Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv",
    "genes": TMP_DIR / 'Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv',
    "obs": TMP_DIR / 'Wu_etal_2021_BRCA_scRNASeq/metadata.csv',
}


# Download the data
print("Downloading data", flush=True)
urllib.request.urlretrieve(
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz', 
    FILE_PATHS["tar_gz"]
)

# Extract the data
scrna_file = tarfile.open(FILE_PATHS["tar_gz"])
scrna_file.extractall(TMP_DIR)
scrna_file.close()

# Read data
print("Reading count matrix", flush=True)
counts = sc.read_mtx(FILE_PATHS["counts_mtx"], dtype='float32')

print("Reading barcodes", flush=True)
barcodes = pd.read_csv(FILE_PATHS["barcodes"], header = None)
barcodes.columns = ['Barcode']
barcodes.set_index('Barcode', inplace=True)
barcodes.rename_axis(None, axis=0, inplace=True)

print("Reading var names", flush=True)
genes = pd.read_csv(FILE_PATHS["genes"], header = None)
genes.columns = ['Gene']
genes.set_index('Gene', inplace=True)
genes.rename_axis(None, axis=0, inplace=True)

print("Reading obs", flush=True)
obs = pd.read_csv(FILE_PATHS["obs"])
obs = obs.iloc[0::,1::]

print("Setting barcodes as obs index", flush=True)
obs.index  = barcodes.index

# Create adata
print("Creating adata", flush=True)
adata = ad.AnnData(X = counts.X.T, var = genes, obs = obs )

adata.layers["counts"] = adata.X
del adata.X

# Rename fields
rename_obs_keys = { 
    "cell_type": "celltype_major", 
    "cell_type_level2": "celltype_minor",
    "cell_type_level3": "celltype_subset",
    "donor_id": "orig.ident",
    "cancer_subtype": "subtype", # TODO: this is currently not in the config yaml - how to handle this?
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
    "dataset_id": "2021Wu_human_breast_cancer_sc", #"GSE176078",
    "assay": "Chromium Single-Cell v2 3’ and 5’ Chemistry Library", # from metadata from GEO GSE176078 #TODO: ontology
    "sex": "female", #TODO: double check
    "tissue": "breast",
    "disease": "breast cancer",
    # "disease_ontology_term_id": "PATO:0000461", #TODO: ontology
    "organism": "Homo sapiens",
    # "organism_ontology_term_id": "NCBITaxon:10090", #TODO: ontology
    "tissue_general": "breast",
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
adata.var["feature_name"] = adata.var_names
# TODO: can we also get ensembl ids? (adata.var["feature_id"]) 

# Uns
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

# Filter for cancer subtypes
adata = adata[adata.obs["cancer_subtype"].isin(par["cancer_subtypes"])]

# Delete files if requested
if not par["keep_files"]:
    print("Removing files", flush=True)
    for file in FILE_PATHS.values():
        if file.exists():
            print("\t...", file, flush=True)
            file.unlink()

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
