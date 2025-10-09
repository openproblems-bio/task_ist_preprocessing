import spatialdata as sd
from tifffile import imwrite
import sys
import numpy as np
import pandas as pd

input_path = sys.argv[1]
output_path_tif = sys.argv[3]
output_path_tsv = sys.argv[2]

sdata = sd.read_zarr(input_path)

transcripts_coord_systems = sd.transformations.get_transformation(sdata['transcripts'], get_all=True).keys()

image_coord_systems = sd.transformations.get_transformation(sdata["morphology_mip"], get_all=True).keys()

print('Transforming transcripts coordinates', flush=True)
transcripts = sd.transform(sdata['transcripts'], to_coordinate_system='global')

image = sdata['morphology_mip']['scale0'].image.compute().to_numpy()
transformation = sdata['morphology_mip']['scale0'].image.transform.copy()

mip2d = image[0] 
imwrite(
    output_path_tif,           # .tif or .ome.tif
    mip2d,
    dtype=mip2d.dtype,      # keeps original bit-depth, e.g. uint8
    compression="zlib")


transcripts_df = transcripts.compute().loc[:,['x', 'y', 'feature_name']]
transcripts_df['x_int'] = transcripts_df.x.astype(int)
transcripts_df['y_int'] = transcripts_df.y.astype(int)

print("Counting transcripts per (gene,x,y)")
agg = (
    transcripts_df.groupby(["feature_name", "x_int", "y_int"], sort=False, observed=True)
      .size()
      .reset_index(name="MIDCounts")
      .rename(columns={"feature_name": "geneID", "x_int": "x", "y_int": "y"})
)
agg['x'] = agg.x - np.min(agg.x)
agg['y'] = agg.y - np.min(agg.y)

agg[agg.MIDCounts!=0].to_csv(output_path_tsv, sep = "\t", index=False)