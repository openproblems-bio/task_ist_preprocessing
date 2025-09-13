import sys
import numpy as np
import anndata as ad

## VIASH START
par = {
    'input_scrnaseq_reference': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
    'output': 'output.h5ad',
    'output_transcript_assignments': 'output_transcript_assignments.zarr',
    'output_qc_col': 'output_qc_col.h5ad',
    'seed': 0,
}
meta = {
    "resources_dir": "src/control_methods/"
}
## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from util import add_layers_obs_var_to_scrnaseq_ref, create_dummy_transcript_assignment_table

# Generate control adata
np.random.seed(par['seed'])

print('Read input_scrnaseq_reference', flush=True)
adata = ad.read_h5ad(par['input_scrnaseq_reference'])

print("Randomise ct annotations", flush=True)
ct_annotations = np.random.permutation(adata.obs["cell_type"])
adata.obs["cell_type"] = ct_annotations

# Generate expected output of dummy processed spatial data
print("Add required layers, obs and var columns for spatial data", flush=True)
add_layers_obs_var_to_scrnaseq_ref(adata)

print("Delete obsm, obsp and varm", flush=True)
del adata.obsm
del adata.varm
del adata.obsp

print("Create dummy transcript assignment table", flush=True)
sdata_transcripts_only = create_dummy_transcript_assignment_table(adata)

print("Create dummy qc column", flush=True)
adata_qc = ad.AnnData(obs=adata.obs[[]])
adata_qc.obs["passed_QC"] = True

# Write outputs
print("Write h5ad", flush=True)
adata.write_h5ad(par['output'], compression='gzip')

print("Write transcripts zarr", flush=True)
sdata_transcripts_only.write(par['output_transcript_assignments'])

print("Write qc column", flush=True)
adata_qc.write_h5ad(par['output_qc_col'], compression='gzip')