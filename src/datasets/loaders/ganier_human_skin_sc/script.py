from pathlib import Path 
import pandas as pd
import anndata as ad
import urllib.request


## VIASH START
par = {
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2024Ganier_human_skin_sc.h5ad",
    "dataset_id": "2024Ganier_human_skin_sc",
    "dataset_name": "2024Ganier_human_skin_sc",
    "dataset_url": "https://spatial-skin-atlas.cellgeni.sanger.ac.uk",
    "dataset_reference": "https://doi.org/10.1073/pnas.2313326120", #TODO: bibtex not doi, also adjust config.vsh.yaml
    "dataset_summary": "This dataset contains scRNA-seq data from healthy human skin.",
    "dataset_description": "This dataset contains scRNA-seq data from healthy human skin.",
    "dataset_organism": "Homo sapiens"
}
meta = {
    "temp_dir": "./temp/datasets/2024Ganier_human_skin_sc",
}
## VIASH END

# helper variables
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_URL = "https://cellgeni.cog.sanger.ac.uk/spatial-skin-atlas/scrnaseq/bcc_and_normal-CG_portal_fat.h5ad"
# cell x gene alternative: https://datasets.cellxgene.cziscience.com/995fa8f9-8ed4-46d3-8d75-115fa06cc787.h5ad
FILE_PATH = TMP_DIR / "bcc_and_normal-CG_portal_fat.h5ad"


# Download the data
print("Downloading data (~1GB)", flush=True)
urllib.request.urlretrieve(FILE_URL, FILE_PATH)

# Read data 
print("Reading data", flush=True)
adata = ad.read_h5ad(FILE_PATH)

# Filter out bcc samples
adata = adata[~adata.obs["01_sample"].str.startswith("bcc")]

# NOTE: only log-normalized counts are available, in the workflow we'll skip the log-normalization step
# This is not optimal, layers "counts" should be the raw counts
adata.layers["counts"] = adata.X
adata.layers["normalized"] = adata.X
del adata.X

# Rename fields
rename_obs_keys = { 
    "cell_type": "04_celltypes", 
    "cell_type_level2": "05_subcelltypes",
    "batch": "01_sample",
    #TODO other fields, note that the cell x gene download source might be the better alternative
}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})

# Add additional information to obs
#TODO: Finish up the terms according to the ontology
#Ontology schema currently (13.03.2025) used in openproblems (CELLxGENE schema v4.0.0): 
#https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md
#(mentioned here: https://openproblems.bio/documentation/reference/openproblems/src-datasets#file-format:-raw-dataset)
store_info = { 
    "dataset_id": "2024Ganier_human_skin_sc", 
    #TODO other fields, note that the cell x gene download source might be the better alternative
    #"assay": "Chromium Single-Cell v2 3’ and 5’ Chemistry Library", # from metadata from GEO GSE176078 #TODO: ontology
    #"sex": "female", #TODO: double check
    #"tissue": "breast",
    #"disease": "breast cancer",
    ## "disease_ontology_term_id": "PATO:0000461", #TODO: ontology
    #"organism": "Homo sapiens",
    ## "organism_ontology_term_id": "NCBITaxon:10090", #TODO: ontology
    #"tissue_general": "breast",
    ## "tissue_general_ontology_term_id": "UBERON:0000955", #TODO: ontology
    #"development_stage": "adult", 
    ## "development_stage_ontology_term_id": "MmusDv:0000110" #TODO: ontology
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
    if FILE_PATH.exists():
        print("\t...", FILE_PATH, flush=True)
        FILE_PATH.unlink()

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
