# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard

import spatialdata as sd
from spatialdata_io import xenium
import shutil
import os

## VIASH START
par = {
    "input": [
        "resources/datasets_raw/10x_fresh_frozen_mouse_brain_replicates/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs",
        "resources/datasets_raw/10x_fresh_frozen_mouse_brain_replicates/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs",
        "resources/datasets_raw/10x_fresh_frozen_mouse_brain_replicates/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    ],
    "replicate_id": [
        "rep1",
        "rep2",
        "rep3",
    ],
    "output": "resources/datasets/10x_xenium/10x_fresh_frozen_mouse_brain_replicates/dataset.zarr"
}
## VIASH END

sdatas = []

for i, input in enumerate(par["input"]):
    print("parsing the data... ", end="", flush=True)
    replicate_id = par["replicate_id"][i]
    sdata = xenium(
        path=input,
        n_jobs=8,
        cell_boundaries=True,
        nucleus_boundaries=True,
        morphology_focus=True,
        cells_as_circles=True,
    )

    # SpatialData object
    # ├── Images
    # │     ├── 'morphology_focus': DataTree[cyx] (1, 33131, 48358), (1, 16565, 24179), (1, 8282, 12089), (1, 4141, 6044), (1, 2070, 3022)
    # │     └── 'morphology_mip': DataTree[cyx] (1, 33131, 48358), (1, 16565, 24179), (1, 8282, 12089), (1, 4141, 6044), (1, 2070, 3022)
    # ├── Labels
    # │     ├── 'cell_labels': DataTree[yx] (33131, 48358), (16565, 24179), (8282, 12089), (4141, 6044), (2070, 3022)
    # │     └── 'nucleus_labels': DataTree[yx] (33131, 48358), (16565, 24179), (8282, 12089), (4141, 6044), (2070, 3022)
    # ├── Points
    # │     └── 'transcripts': DataFrame with shape: (<Delayed>, 8) (3D points)
    # ├── Shapes
    # │     ├── 'cell_boundaries': GeoDataFrame shape: (162033, 1) (2D shapes)
    # │     ├── 'cell_circles': GeoDataFrame shape: (162033, 2) (2D shapes)
    # │     └── 'nucleus_boundaries': GeoDataFrame shape: (162033, 1) (2D shapes)
    # └── Tables
    #       └── 'table': AnnData (162033, 248)
    # with coordinate systems:
    #     ▸ 'global', with elements:
    #         morphology_focus (Images), morphology_mip (Images), cell_labels (Labels), nucleus_labels (Labels), transcripts (Points), cell_boundaries (Shapes), cell_circles (Shapes), nucleus_boundaries (Shapes)


    sdata2_kwargs = {}

    # process images
    sdata2_kwargs["images"] = {
        # morphology_focus or morphology_mip?
        replicate_id + "_image": sdata.images["morphology_focus"],
    }
    sdata2_kwargs["labels"] = {
        replicate_id + "_cell": sdata.labels["cell_labels"],
        replicate_id + "_nucleus": sdata.labels["nucleus_labels"],
    }
    sdata2_kwargs["points"] = {
        replicate_id + "_transcripts": sdata.points["transcripts"],
    }
    sdata2_kwargs["shapes"] = {
        replicate_id + "_cell_boundaries": sdata.shapes["cell_boundaries"],
        replicate_id + "_nucleus_boundaries": sdata.shapes["nucleus_boundaries"],
    }
    sdata2_kwargs["tables"] = {
        replicate_id + "_cell_table": sdata.tables["table"],
    }

    sdata2 = sd.SpatialData(**sdata2_kwargs)
    

    # TODO: 

    sdatas.append(sdata2)


sdata = sd.concatenate(sdatas)

print(sdata)


print("writing the data... ", end="", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sdata.write(par["output"])
