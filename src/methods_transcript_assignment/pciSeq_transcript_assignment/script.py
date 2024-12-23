import numpy as np
import dask
import spatialdata as sd
import txsim as tx
import anndata as ad
import os
import shutil

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_ist': 'resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr',
  'input_segmentation': 'resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr',
  'transcripts_key': 'transcripts',
  'coordinate_system': 'global',
  'output': '../pciSeq_assigned_transcripts.zarr',

  'input_scrnaseq': 'resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad',
  'sc_cell_type_key': 'cell_type',

  'exclude_genes': None,
  'max_iter': 1000,
  'CellCallTolerance': 0.02,
  'rGene': 20,
  'Inefficiency': 0.2,
  'InsideCellBonus': 2,
  'MisreadDensity': 0.00001,
  'SpotReg': 0.1,
  'nNeighbors': 3,
  'rSpot': 2,
  'save_data': False,
  'dtype': np.float64
}
meta = {
  'name': 'pciSeq_transcript_assignment'
}
## VIASH END

# Read input
print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])

# Check if coordinate system is available in input data
transcripts_coord_systems = sd.transformations.get_transformation(sdata[par["transcripts_key"]], get_all=True).keys()
assert par['coordinate_system'] in transcripts_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."
segmentation_coord_systems = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True).keys()
assert par['coordinate_system'] in segmentation_coord_systems, f"Coordinate system '{par['coordinate_system']}' not found in input data."

# Transform transcript coordinates to the coordinate system
print('Transforming transcripts coordinates', flush=True)
transcripts = sd.transform(sdata[par['transcripts_key']], to_coordinate_system=par['coordinate_system'])

# In case of a translation transformation of the segmentation (e.g. crop of the data), we need to adjust the transcript coordinates
trans = sd.transformations.get_transformation(sdata_segm["segmentation"], get_all=True)[par['coordinate_system']].inverse()
transcripts = sd.transform(transcripts, trans, par['coordinate_system'])

# Assign cell ids to transcripts
print('Assigning transcripts to cell ids', flush=True)
y_coords = transcripts.y.compute().to_numpy()
x_coords = transcripts.x.compute().to_numpy()

#Added for pciSeq
#TODO this will immediately break when the name of the gene isn't feature_name
transcripts_dataframe = sdata[par['transcripts_key']].compute()[['feature_name']] 
transcripts_dataframe['x'] = x_coords
transcripts_dataframe['y'] = y_coords

#same as before
label_image = sdata_segm["segmentation"]["scale0"].image.to_numpy() #TODO: mabye this line needs generalization (DataTree vs DataArray)

# Grab all the pciSeq parameters
opts_keys = [#'exclude_genes',
  'max_iter',
  'CellCallTolerance',
  'rGene',
  'Inefficiency',
  'InsideCellBonus',
  'MisreadDensity',
  'SpotReg',
  'nNeighbors',
  'rSpot',
  'save_data']

opts = {k: par[k] for k in opts_keys}

input_scrnaseq = ad.read_h5ad(par['input_scrnaseq'])
input_scrnaseq.X = input_scrnaseq.layers['counts']

assignments, cell_types = tx.preprocessing.run_pciSeq(
    transcripts_dataframe,
    label_image,
    input_scrnaseq,
    par['sc_cell_type_key'],
    opts
)

#assign transcript -> cell
cell_id_dask_series = dask.dataframe.from_dask_array(
    dask.array.from_array(
        assignments['cell'].to_numpy(), chunks=tuple(sdata[par['transcripts_key']].map_partitions(len).compute())
    ), 
    index=sdata[par['transcripts_key']].index
)

sdata[par['transcripts_key']]["cell_id"] = cell_id_dask_series 

# create new .obs for cells based on the segmentation output (corresponding with the transcripts 'cell_id')
cell_types['type'] = cell_types['type'].replace({'None':'None_sp'})
cell_types.insert(0, 'cell_id', cell_types.index)
cell_types.rename(columns={'type':'cell_type','prob':'cell_type_prob'}, inplace=True)

assert 0 not in cell_types['cell_id'], "Found '0' in cell_id column of assingment output cell matrix"

output_table = ad.AnnData(
      obs=cell_types[['cell_id','cell_type','cell_type_prob']],
      var=sdata.tables["table"].var[[]]
    )

# TODO: Also take care of the following cases:
# - segmentation 3D, transcripts 3D
# - segmentation 3D, transcripts 2D
# - segmentation 2D, transcripts 3D

# Subset sdata to transcripts with cell ids

print('Subsetting to transcripts cell id and cell type data', flush=True)
sdata_transcripts_only = sd.SpatialData(
  points={
    "transcripts": sdata[par['transcripts_key']]
  },
  tables={
    "table": output_table
  }
)

print('Write transcripts with cell ids and cell types', flush=True)
if os.path.exists(par["output"]):
  shutil.rmtree(par["output"])
sdata_transcripts_only.write(par['output'])


