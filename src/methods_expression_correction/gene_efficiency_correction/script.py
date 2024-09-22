import anndata as ad
import txsim as tx

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_spatial_with_cell_types': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_with_celltypes.h5ad',
  'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
  'celltype_key': 'cell_type',
  'output': 'spatial_corrected.h5ad',
}
meta = {
  'name': 'gene_efficiency_correction',
}
## VIASH END

# Optional parameter check: For this specific correction method the par['input_sc'] is required
assert par['input_scrnaseq_reference'] is not None, 'Single cell input is required for this expr correction method.'
    
# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input_spatial_with_cell_types'])
adata_sc = ad.read_h5ad(par['input_scrnaseq_reference'])
adata_sp.layers["normalized_uncorrected"] = adata_sp.layers["normalized"]
adata_sp_reduced = adata_sp[:,adata_sc.var_names].copy()
obs_sp = adata_sp.obs.copy()

# Apply gene efficiency correction
print('Annotating cell types', flush=True)
adata_sp_reduced = tx.preprocessing.gene_efficiency_correction(
    adata_sp_reduced, adata_sc, layer_key='normalized', ct_key=par['celltype_key']
)

# Concatenate and reorder #TODO with the assumption that spatial and sc have the same genes things simplify
print('Concatenating and reordering', flush=True)
gene_order = adata_sp.var_names.tolist()
gene_mask = ~adata_sp.var_names.isin(adata_sp_reduced.var_names)
adata_sp = ad.concat([adata_sp[:,gene_mask], adata_sp_reduced], axis=1, join="outer") 
adata_sp = adata_sp[:,gene_order]
adata_sp.obs = obs_sp

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])
