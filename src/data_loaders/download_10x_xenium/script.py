import spatialdata as sd
import anndata as ad
from spatialdata_io import xenium
import shutil
import os

## VIASH START
par = {
    "input": [
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs",
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs",
        "temp/datasets/10x_xenium/2023_10x_mouse_brain_xenium/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    ],
    "replicate_id": [
        "rep1",
        "rep2",
        "rep3",
    ],
    "segmentation_id": [
        "cell",
        "nucleus",
    ],
    "output": "output.zarr",
    "dataset_id": "value",
    "dataset_name": "value",
    "dataset_url": "value",
    "dataset_reference": "value",
    "dataset_summary": "value",
    "dataset_description": "value",
    "dataset_organism": "value",
}
## VIASH END

sdatas = []

for i, input in enumerate(par["input"]):
    print("parsing the data... ", end="", flush=True)
    replicate_id = par["replicate_id"][i]

    # read the data
    sdata = xenium(
        path=input,
        n_jobs=8,
        cells_boundaries=True,
        nucleus_boundaries=True,
        morphology_focus=True,
        cells_as_circles=False,
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


    # rename coordinate system
    sdata.rename_coordinate_systems({"global": replicate_id + "_global"})
    
    # rename images
    sdata.images[replicate_id + "_image"] = sdata.images.pop("morphology_mip")

    # remove morphology_focus
    _ = sdata.imsdata = sd.concatenate(sdatas)
    # rename labels
    sdata.labels[replicate_id + "_cell"] = sdata.labels.pop("cell_labels")
    sdata.labels[replicate_id + "_nucleus"] = sdata.labels.pop("nucleus_labels")

    # rename points
    sdata.points[replicate_id + "_transcripts"] = sdata.points.pop("transcripts")

    # rename shapes
    sdata.shapes[replicate_id + "_cell_boundaries"] = sdata.shapes.pop("cell_boundaries")
    sdata.shapes[replicate_id + "_nucleus_boundaries"] = sdata.shapes.pop("nucleus_boundaries")

    # rename tables
    sdata.tables[replicate_id + "_cell_table"] = sdata.tables.pop("table")

    sdatas.append(sdata)

# concatenate the data
sdata = sd.concatenate(sdatas)

sdata.tables["table"] = ad.AnnData(
    uns={
        "dataset_id": par["dataset_id"],
        "dataset_name": par["dataset_name"],
        "dataset_url": par["dataset_url"],
        "dataset_reference": par["dataset_reference"],
        "dataset_summary": par["dataset_summary"],
        "dataset_description": par["dataset_description"],
        "dataset_organism": par["dataset_organism"],
        "variables": {
            "replicate_id": par["replicate_id"],
            "segmentation_id": par["segmentation_id"],
        }
    }
)

print(sdata)

print("writing the data... ", end="", flush=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sdata.write(par["output"])
