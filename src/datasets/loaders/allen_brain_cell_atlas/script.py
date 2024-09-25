# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
import scipy as sp
import anndata as ad
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

## VIASH START
par = {
    "abca_version": "20230630",
    "regions": ["MB", "TF"],
    "sample_n_obs": 5000,
    "sample_obs_weight": "subclass",
    "sample_transform": "sqrt",
    "sample_seed": None,
    "output": "tmp_dataset.h5ad",
}
meta = {
    "temp_dir": "/tmp/allen_brain_cell_atlas",
}
## VIASH END

if par["sample_seed"]:
    np.random.seed(par["sample_seed"])

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

print("Reading obs", flush=True)
obs = pd.read_csv(
    abc_cache.get_metadata_path(directory="WMB-10X", file_name="cell_metadata_with_cluster_annotation"),
    index_col=0,
    dtype=defaultdict(
        lambda: "category",
        cluster_label="int",
        x="float",
        y="float",
        region_of_interest_order="int"
    )
)

print("Filtering obs based on regions", flush=True)
obs = obs[obs["anatomical_division_label"].isin(REGIONS)]

if par["sample_n_obs"] and par["sample_n_obs"] < obs.shape[0]:
    print("Filtering obs based on n_obs", flush=True)
    col = par["sample_obs_weight"]

    if col:
        weights = obs.groupby(col).size()

        if par["sample_transform"] == "sqrt":
            weights = weights.apply(lambda x: np.sqrt(x))
        elif par["sample_transform"] == "log":
            weights = weights.apply(lambda x: np.log(x))

        obs = obs.sample(n=par["sample_n_obs"], weights=obs[col].map(weights))
    else:
        obs = obs.sample(n=par["sample_n_obs"])


# From abc_cache.list_data_files("WMB-10Xv2")
# TODO: potentially also load other chemistries (currently only 10Xv2)

print("Downloading and reading expression matrices", flush=True)
adatas = []
for region in REGIONS:
    try:
        print(f"Downloading h5ad file for region {region}", flush=True)
        adata_path = abc_cache.get_data_path(directory="WMB-10Xv2", file_name=f"WMB-10Xv2-{region}/raw")

        print(f"Reading h5ad for region {region}", flush=True)
        adata = ad.read_h5ad(str(adata_path))

        # filter cells
        adata = adata[adata.obs_names.isin(obs.index)].copy()

        # add region to obs
        adata.obs["region"] = region

        # move counts to layer
        adata.layers["counts"] = adata.X
        del adata.X
        
        # add anndata to list
        adatas.append(adata)
    except Exception as e:
        print(f"Error reading {region}: {e}")

print("Concatenating data", flush=True)
adata = ad.concat(adatas, merge="first")
del adatas

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
