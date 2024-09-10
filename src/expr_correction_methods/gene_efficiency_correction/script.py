import anndata as ad
import txsim as tx

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'spatial_with_celltypes.h5ad',
  'input_sc': 'sc_norm_counts.h5ad',
  'celltype_key': 'celltype',
  'output': 'spatial_corrected.h5ad',
}
meta = {
  'name': 'gene_efficiency_correction',
}
## VIASH END

    
# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input'])
adata_sc = ad.read_h5ad(par['input_sc'])
adata_sp.layers["lognorm_uncorrected"] = adata_sp.layers["lognorm"]
adata_sp_reduced = adata_sp[:,adata_sc.var_names].copy()
obs_sp = adata_sp.obs.copy()

# Annotate cell types
print('Annotating cell types', flush=True)
adata_sp_reduced = tx.preprocessing.gene_efficiency_correction(
    adata_sp_reduced, adata_sc, layer_key='lognorm', ct_key=par['celltype_key']
)

gene_order = adata_sp.var_names.tolist()
gene_mask = ~adata_sp.var_names.isin(adata_sp_reduced.var_names)
adata_sp = ad.concat([adata_sp[:,gene_mask], adata_sp_reduced], axis=1, join="outer") 
adata_sp = adata_sp[:,gene_order]
adata_sp.obs = obs_sp

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])