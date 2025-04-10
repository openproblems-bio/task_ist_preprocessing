import anndata as ad
import txsim as tx
import scvi
import pandas as pd
import scanpy as sc
import scipy
import numpy as np

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_spatial_with_cell_types': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_with_cell_types.h5ad',
  'celltype_key': 'cell_type',
  'output': '../resolvi_spatial_corrected.h5ad',
  'n_hidden': 32,
  'encode_covariates': False,
  'downsample_counts': True
}
meta = {
  'name': 'resolvi_correction',
}
## VIASH END

# NOTE/TODO: for grid search:
# - n_hidden: 32 (default), 64, 128
# - encode_covariates: False(default)/True
# - downsample_counts: True(default)/False

# Optional parameter check: For this specific correction method the par['input_sc'] is required
   
# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input_spatial_with_cell_types'])
adata_sp.layers["normalized_uncorrected"] = adata_sp.layers["normalized"]

print("Filter cells with <5 counts")
sc.pp.filter_cells(adata_sp, min_genes=5)

spatial_array = np.stack([adata_sp.obs['centroid_x'].values, adata_sp.obs['centroid_y'].values], axis=1)
adata_sp.obsm['X_spatial'] = spatial_array

# Apply gene efficiency correction
print('Running ResolVI', flush=True)

scvi.external.RESOLVI.setup_anndata(adata_sp, labels_key=par['celltype_key'], layer="counts")

supervised_resolvi = scvi.external.RESOLVI(adata_sp, semisupervised=True, 
  n_hidden = par['n_hidden'], 
  encode_covariates = par['encode_covariates'], 
  downsample_counts = par['downsample_counts'])
supervised_resolvi.train(max_epochs=100)

samples_corr = supervised_resolvi.sample_posterior(
        model=supervised_resolvi.module.model_corrected,
        return_sites=['px_rate'],
        summary_fun={"post_sample_q50": np.median},
        num_samples=20, return_samples=False, batch_size=4000) #batch_steps was not a parameter
samples_corr = pd.DataFrame(samples_corr).T

samples = supervised_resolvi.sample_posterior(
    model=supervised_resolvi.module.model_residuals,
    return_sites=[
        'mixture_proportions', 'mean_poisson', 'per_gene_background', 
        'diffusion_mixture_proportion', 'per_neighbor_diffusion', 'px_r_inv'
        ],
    num_samples=20, return_samples=False, batch_size=4000)
samples = pd.DataFrame(samples).T


adata_sp.obsm["X_resolVI"] = supervised_resolvi.get_latent_representation()

# TODO these 2 lines threw errors because 'obs' was not generated in samples_corr
# adata_sp.layers["generated_expression"] = scipy.sparse.csr_matrix(samples_corr.loc['post_sample_q25', 'obs'])
# adata_sp.layers["generated_expression_mean"] = scipy.sparse.csr_matrix(samples_corr.loc['post_sample_means', 'obs'])

adata_sp.layers["corrected_counts"] = adata_sp.layers['counts'].multiply((samples_corr.loc['post_sample_q50', 'px_rate'] / (
    1.0 + samples_corr.loc['post_sample_q50', 'px_rate'] + samples.loc['post_sample_means', 'mean_poisson']))).tocsr()

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])
