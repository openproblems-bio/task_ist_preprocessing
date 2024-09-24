# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

from pathlib import Path
import pandas as pd
import anndata as ad
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import scipy as sp
import scanpy as sc

## VIASH START
par = {
    "abca_version": "20230630",
    "regions": ["CTXsp", "HPF", "HY", "Isocortex-1", "Isocortex-2", "Isocortex-3", "Isocortex-4", "MB", "OLF", "TF"],
    "output": "tmp_dataset.h5ad",
}
meta = {
    "temp_dir": "/tmp/allen_brain_cell_atlas",
}
## VIASH END

# helper variables
VERSION = par["abca_version"]
REGIONS = par["regions"]
TMP_DIR = Path(meta["temp_dir"] or "/tmp")

print("Loading manifest", flush=True)
# saved to TMPDIR / releases/{VERSION}/manifest.json
abc_cache = AbcProjectCache.from_cache_dir(TMP_DIR)
abc_cache.load_manifest(
    f"releases/{VERSION}/manifest.json"
)

print("Downloading metadata", flush=True)
metadata_files = [
    "cell_metadata_with_cluster_annotation",
]
for file in metadata_files:
    abc_cache.get_metadata_path(directory="WMB-10X", file_name=file)

print("Reading obs", flush=True)
obs = pd.read_csv(
    TMP_DIR / f"metadata/WMB-10X/{VERSION}/views/cell_metadata_with_cluster_annotation.csv",
    index_col=0,
)

print("Downloading expression matrices", flush=True)
# From abc_cache.list_data_files("WMB-10Xv2") # TODO: potentially also load other chemistries (currently only 10Xv2)
for region in REGIONS:
    print(f"Downloading h5ad file for region {region}", flush=True)
    file = f"WMB-10Xv2-{region}/raw"
    abc_cache.get_data_path(directory="WMB-10Xv2", file_name=file)

print("Reading expression matrices", flush=True)
adatas = []
for region in REGIONS:
    print(f"Reading h5ad for region {region}", flush=True)
    adata = ad.read_h5ad(
        TMP_DIR / f"expression_matrices/WMB-10Xv2/{VERSION}/WMB-10Xv2-{region}-raw.h5ad"
    )
    adata = adata[adata.obs_names.isin(obs.index)]
    adata.obs["region"] = region
    counts = adata.X
    del adata.X

    # make sure counts is sparse
    if not isinstance(counts, sp.sparse.csr_matrix):
        counts = sp.sparse.csr_matrix(counts)
    adata.layers["counts"] = counts
    
    # add anndata to list
    adatas.append(adata)

print("Concatenating data", flush=True)
adata = ad.concat(adatas, merge="first")
del adatas

print("Filtering data", flush=True)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_genes=1)

print("Processing .obs")
adata.obs = obs.loc[adata.obs.index]

# rename fields
rename_obs_keys = {
    "dataset_id": "dataset_label",
    "assay": "library_method",
    "cell_type": "class", 
    "cell_type_level2": "subclass",
    "cell_type_level3": "supertype",
    "cell_type_level4": "cluster",
    "donor_id": "donor_label",
    "sex": "donor_sex",
    "tissue": "region_of_interest_acronym",
    "batch": "library_label",
    # #TODO "cell_type_unified" (?), maybe call the unified one "cell_type" and the original one "cell_type_level1"
    # other keys: "assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id"
    # "diseases_ontology_term_id", "is_primary_data", "organism_ontology_term_id", "self_reported_ethnicity", 
    # "self_reported_ethnicity_ontology_term_id", "sex_ontology_term_id", "suspension_type", 
    # "suspension_type_ontology_term_id", "tissue_ontology_term_id", "tissue_general_ontology_term_id", "soma_joinid"
}
adata.obs = adata.obs.rename(columns={old:new for new,old in rename_obs_keys.items()})

# add additional information to obs
store_info = {
    "disease": "normal",
    "disease_ontology_term_id": "PATO:0000461",
    "organism": "Mus musculus",
    "organism_ontology_term_id": "NCBITaxon:10090",
    "tissue_general": "brain",
    "tissue_general_ontology_term_id": "UBERON:0000955",
    "development_stage": "adult", # from metadata at GEO GSE246717: all ages >= 51 days
    "development_stage_ontology_term_id": "MmusDv:0000110"
}
for key, value in store_info.items():
    adata.obs[key] = pd.Categorical(value, categories=[value])

# remove undesired columns
for key in adata.obs.columns:
    if (key not in rename_obs_keys.keys()) and (key not in store_info.keys()):
        print(f"Removing .obs['{key}']")
        del adata.obs[key]

# Var
adata.var["feature_id"] = adata.var_names
adata.var = adata.var.rename(columns={"gene_symbol": "feature_name"})

# Uns
for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
    adata.uns[key] = par[key]

print(f"Output: {adata}")

# Write data
adata.write_h5ad(par["output"])
