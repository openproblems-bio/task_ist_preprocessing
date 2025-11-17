from pathlib import Path 
import pandas as pd
import anndata as ad
import scanpy as sc
import urllib.request

## VIASH START
par = {
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2020Travaglini_human_lung_sc.h5ad",
    "dataset_id": "2020Travaglini_human_lung_sc",
    "dataset_name": "2020Travaglini_human_lung_sc",
    "dataset_url": "https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293",
    "dataset_reference": "https://doi.org/10.1038/s41586-020-2922-4", #TODO: bibtex not doi, also adjust config.vsh.yaml
    "dataset_summary": "This dataset contains scRNA-seq data from human lung cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human lung cells.",
    "dataset_organism": "Homo sapiens"
}
meta = {
    "temp_dir": "./temp/datasets/2020Travaglini_human_lung_sc",
}
## VIASH END

# helper variables
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_URL = "https://datasets.cellxgene.cziscience.com/060e8716-9f0e-4773-9417-582ddc9ba7ab.h5ad"
FILE_PATH = TMP_DIR / "060e8716-9f0e-4773-9417-582ddc9ba7ab.h5ad"


# Download the data
print("Downloading data (~5.5GB)", flush=True)
urllib.request.urlretrieve(FILE_URL, FILE_PATH)

adata = ad.read_h5ad(FILE_PATH)

# Subset to Travaglini/Krasnow dataset
adata = adata[adata.obs['dataset'] == "Krasnow_2020"]

# Filter out cell types with less than 30 cells
ct_value_count = adata.obs['ann_finest_level'].value_counts()
cts_to_keep = ct_value_count[ct_value_count > 30].index.tolist()
adata = adata[adata.obs['ann_finest_level'].isin(cts_to_keep)]


adata.layers["counts"] = adata.raw.X
sc.pp.filter_genes(adata, min_counts=1)
del adata.X
del adata.raw

# Remove cell_type obs column as we'll assign ann_finest_level as cell_type
del adata.obs["cell_type"]

# Rename fields
rename_obs_keys = { 
    "cell_type": "ann_finest_level", 
    "batch": "sample",
    "assay": "assay",
    "assay_ontology_term_id": "assay_ontology_term_id",
    "cell_type_ontology_term_id": "cell_type_ontology_term_id",
    "development_stage": "development_stage",
    "development_stage_ontology_term_id": "development_stage_ontology_term_id",
    "disease": "disease",
    "disease_ontology_term_id": "disease_ontology_term_id",
    "is_primary_data": "is_primary_data",
    "self_reported_ethnicity": "self_reported_ethnicity",
    "self_reported_ethnicity_ontology_term_id": "self_reported_ethnicity_ontology_term_id",
    "sex": "sex",
    "sex_ontology_term_id": "sex_ontology_term_id",
    "suspension_type": "suspension_type",
    "tissue": "tissue",
    "tissue_ontology_term_id": "tissue_ontology_term_id",
}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})

# Add additional information to obs
#TODO: Finish up the terms according to the ontology
#Ontology schema currently (13.03.2025) used in openproblems (CELLxGENE schema v4.0.0): 
#https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md
#(mentioned here: https://openproblems.bio/documentation/reference/openproblems/src-datasets#file-format:-raw-dataset)
store_info = { 
    "dataset_id": "2020Travaglini_human_lung_sc", 
    "organism": "Homo sapiens",
    # "organism_ontology_term_id": "NCBITaxon:10090", #TODO: ontology
}
for key, value in store_info.items():
    adata.obs[key] = pd.Categorical([value] * adata.n_obs, categories=[value])

# Remove undesired columns
for key in adata.obs.columns:
    if (key not in rename_obs_keys.keys()) and (key not in store_info.keys()):
        print(f"Removing .obs['{key}']")
        del adata.obs[key]

# Var
adata.var["gene_symbol"] = adata.var["feature_name"]
adata.var["feature_id"] = adata.var_names
adata.var_names = adata.var["feature_name"]
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()
adata.var.index.name = None
adata.var["feature_name"] = adata.var_names.astype(str)

# Uns
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

# Delete files if requested
if not par["keep_files"]:
    print("Removing files", flush=True)
    if FILE_PATH.exists():
        print("\t...", FILE_PATH, flush=True)
        FILE_PATH.unlink()

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
