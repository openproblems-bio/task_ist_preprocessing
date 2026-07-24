import os
import spatialdata as sd
import sopa
import anndata as ad
import pandas as pd
import numpy as np
import dask
import dask.dataframe as dd

## VIASH START
par = {
    "input_ist": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
    "input_segmentation": "resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr",
    "transcripts_key": "transcripts",
    "coordinate_system": "global",
    "output": "temp/comseg/transcripts.zarr",
    
    "patch_width": 2000,
    "patch_overlap": 50,
    "transcript_patch_width": 200,
    "mean_cell_diameter": 15.0,
    "max_cell_radius": 25.0,
    "alpha": 0.5,
    "min_rna_per_cell": 5,
    "gene_column": "feature_name",
    "norm_vector": False,
    "allow_disconnected_polygon": True,
    "n_workers": 4,
    "worker_memory_limit": "0",
}
meta = {
    "name": "comseg",
    "cpus": 15,
}
## VIASH END


# Read input files
print('Reading input files', flush=True)
sdata = sd.read_zarr(par['input_ist'])
sdata_segm = sd.read_zarr(par['input_segmentation'])

# Multi-partition parquet restarts its index at 0 per partition, producing duplicate
# index labels that break sopa's cell_id assignment and the final write with
# "cannot reindex on an axis with duplicate labels". Rebuild the transcripts as a
# single-partition frame with a clean RangeIndex before any sopa op touches them.
_tx = sdata[par['transcripts_key']]
_tx_reset = dd.from_pandas(_tx.compute().reset_index(drop=True), npartitions=1)
_tx_reset.attrs.update(_tx.attrs)
del sdata[par['transcripts_key']]
sdata[par['transcripts_key']] = _tx_reset

# Convert the prior segmentation to polygons
sdata["segmentation_boundaries"] = sd.to_polygons(sdata_segm["segmentation"])
del sdata["segmentation_boundaries"]["label"] # make_transcript_patches will create a new label column and fails if one exists.

# Make patches
sopa.make_image_patches(sdata, image_key="image", patch_width=par["patch_width"], patch_overlap=par["patch_overlap"])

transcript_patch_args = {
    "sdata": sdata,
    "write_cells_centroids": True,
    "patch_width": par["transcript_patch_width"],
    "prior_shapes_key": "segmentation_boundaries",
    "points_key": par["transcripts_key"],
}

sopa.make_transcript_patches(**transcript_patch_args)

# Run ComSeg
config = {
    "dict_scale": {"x": 1, "y": 1, "z": 1},
    "mean_cell_diameter": par["mean_cell_diameter"],
    "max_cell_radius": par["max_cell_radius"],
    "norm_vector": par["norm_vector"],
    "alpha": par["alpha"], 
    "allow_disconnected_polygon": par["allow_disconnected_polygon"],
    "min_rna_per_cell": par["min_rna_per_cell"],
    "gene_column": par["gene_column"],
}


# ComSeg processes each transcript patch independently and is pure-Python, so it
# runs single-threaded unless a parallelization backend is enabled. Use the dask
# backend to spread patches across workers, but CAP concurrency: comseg builds an
# in-memory graph per patch, so peak memory scales with the number of concurrent
# workers. Running one worker per CPU (14 on midcpu) multiplied peak RAM ~14x and
# got the job OOM-killed (exit 137 = cgroup SIGKILL). Bounding the worker count is
# the real memory safeguard (the per-worker memory monitor is neutralised below,
# see the note there).
cpu_cap = max((meta["cpus"] or os.cpu_count() or 1) - 1, 1)
n_workers = max(min(par["n_workers"], cpu_cap), 1)

# `distributed` defaults worker startup to "spawn", which re-imports this script
# in every worker. As a plain viash script it has no `if __name__ == "__main__"`
# guard, so the re-import re-runs the module-level code (creating yet another
# cluster) and Python's spawn bootstrap check aborts with "An attempt has been
# made to start a new process before ... bootstrapping". Force "fork" so workers
# inherit the already-initialised interpreter instead of re-importing it; this is
# the Linux-native default and the mode comseg ran under before the deps bumped.
dask.config.set({"distributed.worker.multiprocessing-method": "fork"})

sopa.settings.parallelization_backend = "dask"
sopa.settings.dask_client_kwargs["n_workers"] = n_workers
sopa.settings.dask_client_kwargs["threads_per_worker"] = 1  # CPU-bound work, avoid GIL contention

# ComSeg's per-patch memory is mostly UNMANAGED (native numba/shapely/graph
# allocations, plus the forked interpreter's copy-on-write pages), which dask
# cannot spill. A per-worker memory_limit therefore never spills here; it only
# lets the nanny kill workers that cross ~95% of the limit -> KilledWorker, even
# though nothing is leaking. dask's "auto" is worse: it divides detected RAM by
# the CPU count (not n_workers), so on an 8-core box each worker gets ~1/8 of RAM
# and a normal ~1 GB patch already trips 95%. We therefore neutralise the monitor
# by default, but sopa reads worker.memory_manager.memory_limit and compares it to
# 4 GiB, so the value must be a real number (0/None makes sopa crash). Setting the
# limit to the FULL detected container memory does both: sopa sees a number, and
# the nanny only ever fires near 95% of the whole cgroup — where the OS + Nextflow
# retry already take over — never prematurely for a single patch. Total memory is
# bounded by n_workers instead. An explicit value (e.g. "20GB") is honoured as-is.
_mem_limit_arg = (par["worker_memory_limit"] or "").strip().lower()
if _mem_limit_arg in ("", "0", "none", "off", "disabled"):
    from distributed.system import MEMORY_LIMIT as _CONTAINER_MEMORY
    worker_memory_limit = _CONTAINER_MEMORY
else:
    worker_memory_limit = par["worker_memory_limit"]
sopa.settings.dask_client_kwargs["memory_limit"] = worker_memory_limit

print(
    f"Running ComSeg with dask backend, n_workers={n_workers} "
    f"(cap {cpu_cap}), memory_limit={worker_memory_limit!r}",
    flush=True,
)

sopa.segmentation.comseg(sdata, config)

# Assign transcripts to cell ids
sopa.spatial.assign_transcript_to_cell(
    sdata,
    points_key="transcripts",
    shapes_key="comseg_boundaries",
    key_added="cell_id",
    unassigned_value=0
)

# Create output SpatialData 

# Create objects for cells table
print('Creating objects for cells table', flush=True)
unique_cells = np.unique(sdata["transcripts"]["cell_id"])
zero_idx = np.where(unique_cells == 0)
if len(zero_idx[0]): 
    unique_cells=np.delete(unique_cells, zero_idx[0][0])

cell_id_col = pd.Series(unique_cells, name='cell_id', index=unique_cells)

# Create transcripts only sdata
print('Subsetting to transcripts cell id data', flush=True)
sdata_transcripts_only = sd.SpatialData(
    points={
        "transcripts": sdata['transcripts']
    },
    tables={
        "table": ad.AnnData(
          obs=pd.DataFrame(cell_id_col),
        )
    }
)


output_path = par['output']
sdata_transcripts_only.write(output_path, overwrite=True)


