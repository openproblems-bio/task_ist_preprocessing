import anndata as ad
import tangram as tg
import torch

## VIASH START
par = {
    'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
    'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
    'output': 'spatial_with_celltypes.h5ad',
    'celltype_key': 'cell_type',
    'mode': 'cells',
    'num_epochs': 1000,
}
meta = {
    'name': 'tangram',
}
## VIASH END

# GPU check
if torch.cuda.is_available():
    device = "cuda:0"
else:
    device = "cpu"

# Optional parameter check: For this specific annotation method the par['input_spatial_normalized_counts'] and par['input_scrnaseq_reference'] are required
assert par['input_spatial_normalized_counts'] is not None, 'Spatial input is required for this annotation method.'
assert par['input_scrnaseq_reference'] is not None, 'Single cell input is required for this annotation method.'

# Read input
adata_sp = ad.read_h5ad(par['input_spatial_normalized_counts'])
adata_sc = ad.read_h5ad(par['input_scrnaseq_reference'])

# use log1p noramlized values  
adata_sc.X = adata_sc.layers['normalized']
adata_sp.X = adata_sp.layers['normalized']
    
adata_sp_orig = adata_sp.copy()

# use all the genes from adata_sp as markers for tangram
markers = adata_sp.var_names.tolist()
    
# Removes genes that all entries are zero. Finds the intersection between adata_sc, adata_st and given marker gene list, 
# save the intersected markers in two adatas. Calculates density priors and save it with adata_st
tg.pp_adatas(adata_sc=adata_sc, adata_sp=adata_sp, genes=markers)
    
# Map single cell data (`adata_sc`) on spatial data (`adata_sp`).
# density_prior (str, ndarray or None): Spatial density of spots, when is a string, value can be 'rna_count_based' or 
# 'uniform', when is a ndarray, shape = (number_spots,). 
# use 'uniform' if the spatial voxels are at single cell resolution (e.g. MERFISH). 'rna_count_based', assumes that 
# cell density is proportional to the number of RNA molecules.
adata_map = tg.map_cells_to_space(
    adata_sc=adata_sc,
    adata_sp=adata_sp,
    device=device,
    mode=par['mode'],
    num_epochs=par['num_epochs'],
    density_prior='uniform'
)
    
# Spatial prediction dataframe is saved in `obsm` `tangram_ct_pred` of the spatial AnnData
tg.project_cell_annotations(
    adata_map = adata_map,
    adata_sp = adata_sp, 
    annotation=par['celltype_key']
)

# Use original without extra layers generated from tangram
df = adata_sp.obsm['tangram_ct_pred'].copy()
adata_sp = adata_sp_orig.copy()

# Set the cell type annotation
adata_sp.obs[par['celltype_key']] = df.idxmax(axis=1)


# # Normalize by row before setting the score
# normalized_df = df.div(df.sum(axis=1), axis=0)
# max_values = normalized_df.max(axis=1)
# adata_sp.obs['tangram_score'] = max_values
# adata_sp.obsm['ct_tangram_scores'] = normalized_df

# Write output
adata_sp.write_h5ad(par['output'])