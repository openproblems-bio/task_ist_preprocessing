import pandas as pd
import anndata as ad
import spatialdata as sd
import txsim as tx

## VIASH START
par = {
  'input': 'resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr',
  'alpha': 0.0,
  'output': 'cell_volumes.h5ad',
}
## VIASH END

print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input'])

print('Determine cell ids', flush=True)
cell_ids = sorted(sdata["transcripts"]["cell_id"].unique())
if cell_ids[0] == 0:
    cell_ids = cell_ids[1:]

print('Init AnnData object', flush=True)
adata = ad.AnnData(
  obs=pd.DataFrame(
    index=cell_ids,
    data={"cell_id": cell_ids}
  ),
  uns={"spots": sdata["transcripts"].compute()[["x","y","cell_id"]]}
)

print('Calculate alpha shape area', flush=True)
tx.preprocessing.calculate_alpha_area(adata=adata, alpha=par['alpha'], cell_id_col="cell_id")
adata.obs["volume"] = adata.obs["alpha_area"]

# TODO: Currently we actually calculate the area instead of the volume. Also we ignore the different z positions/layers.
#       We could improve the calculation with 1. taking into account the scaling in z direction and 2. calculating
#       polygons for each z layer and then summing up the volumes of the polygons.
#       Another option would be an alpha complex instead of an alpha shape. Ideally the x,y,z coordinates are 
#       transformed into physical space then first. There's a 3d option in https://pypi.org/project/alphashape/

print('Write cell volumes h5ad', flush=True)
adata.write_h5ad(par['output'])
