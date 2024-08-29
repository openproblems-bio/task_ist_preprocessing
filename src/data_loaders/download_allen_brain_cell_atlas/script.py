# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

from pathlib import Path
import pandas as pd
import anndata as ad
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# env setup:
# pip install -U git+https://github.com/alleninstitute/abc_atlas_access

## VIASH START
par = {
    "version": "20230630",
    "regions": ["OLF", "TH"],
    "output": f"abc_atlas_20230630.h5ad",
}
meta = {
    "temp_dir": "/tmp",
}
## VIASH END

# helper variables
VERSION = par["version"]
REGIONS = par["regions"]
TMP_DIR = Path(meta["temp_dir"] or "/tmp")

abc_cache = AbcProjectCache.from_cache_dir(TMP_DIR)
abc_cache.load_manifest(
    f"releases/{VERSION}/manifest.json"
)  # saved to TMPDIR / releases/{VERSION}/manifest.json

# From abc_cache.list_data_files('WMB-10Xv2') # TODO: potentially also load other chemistries (currently only 10Xv2)
count_matrix_files = [f'WMB-10Xv2-{region}/raw' for region in REGIONS]

# From abc_cache.list_metadata_files('WMB-10Xv2')
metadata_files = [
    'cell_metadata_with_cluster_annotation',
    #'gene',
    #'region_of_interest_metadata'
]

# Download data
for file in count_matrix_files:
    abc_cache.get_data_path(directory='WMB-10Xv2', file_name=file)

for file in metadata_files:
    abc_cache.get_metadata_path(directory='WMB-10X', file_name=file)

# Read an concatenate the data
obs = pd.read_csv(
    TMP_DIR
    / f"metadata/WMB-10X/{VERSION}/views/cell_metadata_with_cluster_annotation.csv",
    index_col=0,
)

adatas = []
for region in REGIONS:
    adata = ad.read_h5ad(
        TMP_DIR / f"expression_matrices/WMB-10Xv2/{VERSION}/WMB-10Xv2-{region}-raw.h5ad"
    )
    adata = adata[adata.obs_names.isin(obs.index)]
    adata.obs["region"] = region
    adatas.append(adata)

adata = ad.concat(adatas)

# Renaming etc. to match the api

# Layers
adata.layers["counts"] = adata.X

# Obs
new_to_old_obs_keys = {
    "dataset_id": "dataset_label", "assay": "library_method", "cell_type":'class', 
    "cell_type_level2": "subclass", "cell_type_level3": "supertype", "cell_type_level4": "cluster",
    "donor_id":'donor_label', "sex": "donor_sex", "tissue": "region_of_interest_acronym", "batch": "library_label",
    # #TODO "cell_type_unified" (?), maybe call the unified one "cell_type" and the original one "cell_type_level1"
    # other keys: "assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id"
    # "diseases_ontology_term_id", "is_primary_data", "organism_ontology_term_id", "self_reported_ethnicity", 
    # "self_reported_ethnicity_ontology_term_id", "sex_ontology_term_id", "suspension_type", 
    # "suspension_type_ontology_term_id", "tissue_ontology_term_id", "tissue_general_ontology_term_id", "soma_joinid"
 }
new_key_to_value = {
    "disease": "healthy", "organism": "Mus musculus", "tissue_general": "brain", 
    "development_stage": "adult", # from metadata at GEO GSE246717: all ages >= 51 days
}

adata.obs = obs.rename(columns={old:new for new,old in new_to_old_obs_keys.items()})
for key, value in new_key_to_value.items():
    adata.obs[key] = value
for key in adata.obs.columns:
    if (key not in new_to_old_obs_keys.keys()) and (key not in new_key_to_value.keys()):
        del adata.obs[key]

# Var
adata.var["feature_id"] = adata.var_names
adata.var = adata.var.rename(columns={"gene_symbol":"feature_name"})
adata.var_names = adata.var["feature_name"]
adata.var_names_make_unique()
adata.var.index.name = None

# Uns
adata.uns["dataset_id"] = "2023_Yao_mouse_brain_scRNAseq_10Xv2"
adata.uns["dataset_name"] = "2023_Yao_mouse_brain_scRNAseq_10Xv2"
adata.uns["dataset_url"] = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717"
adata.uns["dataset_reference"] = "10.1038/s41586-023-06812-z"
adata.uns["dataset_summary"] = "A high-resolution scRNAseq atlas of cell types in the whole mouse brain"
adata.uns["dataset_description"] = "See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset."
adata.uns["dataset_organism"] = "Mus musculus"

# Write data
adata.write_h5ad(par["output"])

# Delete the temporary files and directories
import shutil
shutil.rmtree(TMP_DIR)
