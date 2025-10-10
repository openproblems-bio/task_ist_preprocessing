import sys
import anndata as ad
import spatialdata as sd

## VIASH START
par = {
    'input': "resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_processed_complete.zarr",
    'output': "quality_metrics.h5ad",
}
meta = {
    'name': 'quality_metrics',
    "resources_dir": "src/metrics/quality"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from quality_metrics import proportion_of_assigned_reads, proportion_of_annotated_cells

print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input'])

print('Compute metrics', flush=True)
metrics = {
    "proportion_of_assigned_reads": proportion_of_assigned_reads(sdata)[0],
    "proportion_of_annotated_cells": proportion_of_annotated_cells(sdata),
}

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'metric_ids': list(metrics.keys()),
    'metric_values': list(metrics.values())
  }
)
output.write_h5ad(par['output'], compression='gzip')