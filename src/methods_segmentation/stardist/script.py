import os
import shutil
from pathlib import Path
import numpy as np
# numpy>=1.24 removed the deprecated scalar-type aliases (np.bool, np.int,
# np.float, np.long, ...), but stardist/csbdeep still reference them. TensorFlow
# 2.17 forces numpy>=1.26, so we can't downgrade numpy far enough to get them
# back — restore the aliases instead. Each mapped to its Python builtin / numpy
# type as numpy itself did before removal.
for _alias, _target in {
    "bool": bool, "int": int, "float": float, "complex": complex,
    "object": object, "str": str, "long": int, "unicode": str,
}.items():
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)
import xarray as xr
import spatialdata as sd
import anndata as ad
#from csbdeep.utils import normalize
from csbdeep.data import Normalizer, normalize_mi_ma
from stardist.models import StarDist2D



def convert_to_lower_dtype(arr):
    max_val = arr.max()
    if max_val <= np.iinfo(np.uint8).max:
        new_dtype = np.uint8
    elif max_val <= np.iinfo(np.uint16).max:
        new_dtype = np.uint16
    elif max_val <= np.iinfo(np.uint32).max:
        new_dtype = np.uint32
    else:
        new_dtype = np.uint64

    return arr.astype(new_dtype)

## VIASH START
par = {
  "input": "resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr",
  "output": "temp/stardist/segmentation.zarr",
  "model": "2D_versatile_fluo",
  "prob_thresh": None,
  "nms_thresh": None,
  "scale": None,
  "max_object_diameter": None,
}

## VIASH END


# Read image and its transformation
sdata = sd.read_zarr(par["input"])
image = sdata['image']['scale0'].image.compute().to_numpy()
transformation = sdata['image']['scale0'].image.transform.copy()

# Segment image

# Load pretrained model
model = StarDist2D.from_pretrained(par['model'])

# Segment on normalized image 
#labels, _ = model.predict_instances(normalize(image)[0,:,:]) # scale = None, **hyperparams)

# from https://github.com/stardist/stardist/blob/main/examples/other2D/predict_big_data.ipynb
class MyNormalizer(Normalizer):
    def __init__(self, mi, ma):
            self.mi, self.ma = mi, ma
    def before(self, x, axes):
        return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)
    def after(*args, **kwargs):
        assert False
    @property
    def do_after(self):
        return False

mi, ma = np.percentile(image, [1,99.8])
normalizer = MyNormalizer(mi, ma)

# Tunable knobs forwarded through predict_instances_big -> predict_instances.
# A value left as None (i.e. omitted from par) means "use the model's own optimized
# value": thresholds.json for prob_thresh/nms_thresh, no rescaling for scale. This
# keeps the default (no-args) call identical to the pre-tuning behaviour.
eval_params = {
    k: par[k]
    for k in ("prob_thresh", "nms_thresh", "scale")
    if par.get(k) is not None
}
print(f"predict_instances_big overrides: {eval_params}", flush=True)

# predict_instances_big tiles the image and stitches the per-block predictions. Its
# stitching invariant is that EVERY predicted object is smaller than `min_overlap`;
# an object spanning a block seam that is larger than the overlap can't be assigned to
# a single block, which raises "Found object of shape (...), which violates the
# assumption of being smaller than 'min_overlap'". So `min_overlap` must be
# OBJECT-SIZE-based, not block-size-based — the old `block_size // 5.5` shrank the
# overlap below real nuclei/blobs on small panels (min_overlap fell to 64 px while
# objects reached ~110 px). Objects are measured in ORIGINAL image pixels
# (predict_instances undoes `scale` internally), so this bound is in original px and is
# independent of `scale`. `context` is only the receptive-field margin discarded around
# each block, so deriving it from the image size is fine.
block_size = min(image.shape[1] // 3, 4096)
context = int(min(block_size // 5.5, 128))
min_overlap = int(par.get("max_object_diameter") or 192)  # px; must exceed largest object
# predict_instances_big asserts: min_overlap + 2*context < block_size.
block_size = max(block_size, min_overlap + 2 * context + 1)

# Self-heal: if a rare oversized blob (merged nuclei / debris) still exceeds
# `min_overlap`, double it (and grow block_size to keep the geometry constraint) and
# retry, rather than failing the whole segmentation.
while True:
    try:
        labels, _ = model.predict_instances_big(
            image[0, :, :], axes='YX', block_size=block_size,
            min_overlap=min_overlap, context=context,
            normalizer=normalizer, **eval_params,  # n_tiles left to block_size
        )
        break
    except RuntimeError as e:
        if "min_overlap" not in str(e) or min_overlap >= 2048:
            raise
        min_overlap *= 2
        block_size = max(block_size, min_overlap + 2 * context + 1)
        print(
            "predict_instances_big: an object exceeded min_overlap; retrying with "
            f"min_overlap={min_overlap}, block_size={block_size}",
            flush=True,
        )



# Create output
sd_output = sd.SpatialData()
labels = convert_to_lower_dtype(labels)
labels_array = xr.DataArray(labels, name=f'segmentation', dims=('y', 'x'))
parsed_labels = sd.models.Labels2DModel.parse(labels_array, transformations=transformation)
sd_output.labels['segmentation'] = parsed_labels

metadata = sdata.tables["metadata"]
# cell_id is required downstream. Standard Xenium exports carry an explicit
# "cell_id" column, but some exports (e.g. the Xenium WTA preview used for the
# Atera dataset) don't — there the per-cell identifier lives in the table's
# instance_key column (falling back to the obs index).
instance_key = metadata.uns.get("spatialdata_attrs", {}).get("instance_key")
if "cell_id" in metadata.obs.columns:
    cell_id = metadata.obs["cell_id"].values
elif instance_key and instance_key in metadata.obs.columns:
    cell_id = metadata.obs[instance_key].values
else:
    cell_id = metadata.obs.index.values
obs = metadata.obs[[]].copy()
obs["cell_id"] = cell_id
if "region" in metadata.obs.columns:
    obs["region"] = metadata.obs["region"].values
sd_output.tables['table'] = ad.AnnData(
      obs=obs,
      var=metadata.var[[]]
    )

print("Writing output", flush=True)
Path(par["output"]).parent.mkdir(parents=True, exist_ok=True)
if os.path.exists(par["output"]):
    shutil.rmtree(par["output"])
sd_output.write(par["output"])

