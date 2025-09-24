import anndata as ad
import spatialdata as sd
import txsim as tx

## VIASH START
par = {
  'input_spatial_normalized_counts': 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad',
  'input_transcript_assignments': 'resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr',
  'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
  'um_per_pixel': 0.5,
  'output': 'spatial_with_celltypes.h5ad',
  'celltype_key': 'celltype',
}
meta = {
  'name': 'ssam',
}
## VIASH END

# Optional parameter check: For this specific annotation method the par['input_transcript_assignments'] and par['input_sc'] are required
assert par['input_transcript_assignments'] is not None, 'Transcripts input is required for this annotation method.'
assert par['input_scrnaseq_reference'] is not None, 'Single cell input is required for this annotation method.'

# Read input
print('Reading input files', flush=True)
adata_sp = ad.read_h5ad(par['input_spatial_normalized_counts'])
transcripts = sd.SpatialData.read(par['input_transcript_assignments'])['transcripts']
adata_sc = ad.read_h5ad(par['input_scrnaseq_reference'])
adata_sc.X = adata_sc.layers["normalized"]
shared_genes = [g for g in adata_sc.var_names if g in adata_sp.var_names]
adata_sp = adata_sp[:,shared_genes] #TODO: do we want to do this earlier? Maybe not, thinking about transcripts related quality metrics

# Annotate cell types
print('Annotating cell types', flush=True)
adata_sp = tx.preprocessing.run_ssam(
    adata_sp, transcripts.compute(), adata_sc, um_p_px=par['um_per_pixel'], 
    cell_id_col='cell_id', gene_col='feature_name', sc_ct_key=par['celltype_key']
)
adata_sp.obs["cell_type"] = adata_sp.obs["ct_ssam"].astype(str)
#TODO: check if coordinate systems of adata_sp, transcripts and the ssam procedure align and what to do with um_p_px
#      currently ssam most likely outputs bad results since the transcripts are provided in physical space instead of pixel space
#TODO: How to handle 'None_sc' this could annotate spatial cell types as None_sc, these should then be 'None_sp' such that their ignored in the metrics.

# Write output
print('Writing output', flush=True)
adata_sp.write(par['output'])