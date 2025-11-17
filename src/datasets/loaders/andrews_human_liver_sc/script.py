from pathlib import Path 
import pandas as pd
import anndata as ad
import scanpy as sc
import urllib.request

## VIASH START
par = {
    "keep_files": True, # wether to delete the intermediate files
    "output": "./temp/datasets/2022Andrews_human_liver_sc.h5ad",
    "dataset_id": "2022Andrews_human_liver_sc",
    "dataset_name": "2022Andrews_human_liver_sc",
    "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185477",
    "dataset_reference": "https://doi.org/10.1002/hep4.1854",
    "dataset_summary": "This dataset contains scRNA-seq data from human liver cells.",
    "dataset_description": "This dataset contains scRNA-seq data from human liver cells.",
    "dataset_organism": "Homo sapiens"
}
meta = {
    "temp_dir": "./temp/datasets/2022Andrews_human_liver_sc",
}
## VIASH END

# helper variables
TMP_DIR = Path(meta["temp_dir"] or "/tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)
FILE_URL = "https://datasets.cellxgene.cziscience.com/9a3c5218-afab-4023-8139-ccd6f049b6a1.h5ad"
FILE_PATH = TMP_DIR / "9a3c5218-afab-4023-8139-ccd6f049b6a1.h5ad"

# Download the data
print("Downloading data", flush=True)
urllib.request.urlretrieve(FILE_URL, FILE_PATH)

# Read data
print("Reading data", flush=True)
adata = ad.read_h5ad(FILE_PATH)

# Filter out normal samples and unknown cell types
adata = adata[(adata.obs["disease"] == "normal") & (adata.obs["cell_type"] != "unknown")]

# Set layer "counts" and remove adata.X
adata.layers["counts"] = adata.raw.X.copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=1)

# Delete .X and .raw
del adata.X
del adata.raw

# Add additional information to obs
store_info = { 
    "dataset_id": "2022Andrews_human_liver_sc", 
}
for key, value in store_info.items():
    adata.obs[key] = pd.Categorical([value] * adata.n_obs, categories=[value])

# Remove undesired columns
# TODO: could remove some columns (we did this for other datasets)

# Var
adata.var["gene_symbol"] = adata.var["feature_name"]
adata.var["ensembl_id"] = adata.var_names
adata.var_names = adata.var["feature_name"]

# Uns
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

# Delete files if requested
if not par["keep_files"]:
    print("Removing files", flush=True)
    if Path(FILE_PATH).exists():
        print("\t...", FILE_PATH, flush=True)
        Path(FILE_PATH).unlink()

print("Writing adata", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
