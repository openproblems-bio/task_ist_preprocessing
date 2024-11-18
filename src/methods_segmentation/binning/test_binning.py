import subprocess
import os
import numpy as np
import yaml
import spatialdata as sd
import anndata as ad
import shutil
import numpy as np
from spatialdata.models import Labels2DModel
import xarray as xr
import datatree as dt

## VIASH START
meta = {
    "executable": "target/executable/binning/binning",
    "config": "src/methods_segmentation/binning/config.vsh.yaml",
    "resources_dir": "resources"
}
## VIASH END
input_path =  "foo.txt" 
print(">>> Create input test file")
with open(input_path, "w") as file:
  file.write(content)


def run_component(cmd):
    print(f">> Running script as test", flush=True)
    out = subprocess.run(cmd)

    assert out.returncode == 0, f"Script exited with an error. Return code: {out.returncode}"   

def check_input_files(arguments):
    print(">> Checking whether input files exist", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "input" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Input file '{arg['value']}' does not exist"

def check_output_files(arguments):
    print(">> Checking whether output file exists", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "output" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Output file '{arg['value']}' does not exist"

    print(">> Reading h5ad files and checking formats", flush=True)
    for arg in arguments:
        if arg["type"] != "file" or arg["direction"] != "output":
            continue
        check_format(arg)



print(">>> Test finished successfully")