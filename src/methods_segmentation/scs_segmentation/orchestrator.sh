#!/bin/bash

## VIASH START
par_input= "../task_ist_preprocessing/resources_test/common/2023_10x_mouse_brain_xenium/dataset.zarr"
par_output= "segmentation.zarr"
meta_resources_dir="../task_ist_preprocessing/src/methods_segmentation/scs"
## VIASH END

par_intermediate_dir=$(mktemp -d -p "$(pwd)" tmp-processing-XXXXXXXX)
echo "meta dir: $meta_resources_dir"

# Access the YAML file
CONDA_ENV_FILE="$meta_resources_dir/environment.yaml"


echo "running SCS orchestrator"

# Create intermediate directory
mkdir -p "$par_intermediate_dir"

which python
# Step 1: Run Python script to reformat input in the first Python environment
 python "input.py" \
  "$par_input" "$par_intermediate_dir/temp.tsv" "$par_intermediate_dir/temp.tif"

export CONDA_DIR=/opt/conda
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
# Put conda in path so we can use conda activate
export PATH=$CONDA_DIR/bin:$PATH

# Create and activate the second Python environment
# Initialize conda for bash
eval "$(/opt/conda/bin/conda shell.bash hook)"
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda env create -f $CONDA_ENV_FILE

conda activate scs_tf

pip install spateo-release
pip install scanpy==1.9.3
pip install tensorflow_addons==0.18.0
pip install --user scikit-misc==0.3.1
pip install --user numpy==1.23.5

#clone SCS repo
git clone https://github.com/chenhcs/SCS.git 


mv script.py SCS/
cd SCS

# Step 2: Run SCS in the second Python environment
python "script.py" \
  "$par_intermediate_dir/temp.tsv" "$par_intermediate_dir/temp.tif" 

conda deactivate

which python ## testing the versions are ok again

cd ../
# Step 3: Run output reformatting in the first Python environment
/usr/local/bin/python "output.py" \
"$par_input" "SCS/results/spot2cell_0:0:0:0.txt" "$par_output"
