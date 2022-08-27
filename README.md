# txsim-pipeline
Pipeline using txsim to compare single cell and spatial transcriptomics data


# Installation
In a clean `conda` environment, run in the terminal:

```git clone https://github.com/theislab/txsim-pipeline.git```

In your `conda` environment, use the `conda-forge` channel to install:

- `mamba`
- `snakemake`

Note: You can also use the `bioconda` channel to install snakemake

## Dependencies

This pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) which will generate the conda environments needed for each job. However, the `txsim` package is still in developement at this time. In order to properly use the pipeline, follow these steps:
- Clone the [txsim repo](https://github.com/theislab/txsim.git) (Note: do not clone inside the `txsim-pipeline` folder) 
- In each folder in `envs`, change the file path `"/mnt/c/..."` to the path of the cloned `txsim` repo


# Running the pipeline

## Setup
Before you can run your pipeline, make sure to edit the config file

Most important is to change all the paths to your desired `test` data and results folder.

Once this is done, try running:

```snakemake -n``` 

for a dry run to ensure everything was installed correctly. Since this only calculates the dependencies and does not actually run any jobs, you can usually run this commnad on a submit node for a cluster if so desired.

## Adding data
To add or rename a dataset, simply follow the convention of `test` for the folder, image, spots, and scRNAseq data paths. The `segmented_image` path is optional if one has manually segmented the nuclei (i.e. using ImageJ).

## Changing the workflow
All supported workflows are shown in the `full_config.yaml` file in `configs`. The names of the batches are arbitrary. The names of the method groups (such as `segmentation`) are also arbitrary but must match the names in the `workflow` key. For each method group, group-wide parameters should be put in the key `[GROUP]_params`. The names of the methods are not arbitrary and are case sensitive. In general, the format of new parameters should follow that of the `full_config.yaml`.

All parameters should be entered as a list, even if it is 1 value. If a method takes a dictionary as a parameter, the config should have a list of dictionaries. This can also be done using the "-" character for a multi-line list

Note: the `custom` method in segmentation will use the `segmented_image` key in `DATA_SCENARIOS`

## Details on combinations

The ConfigParser will generate all possible combinations of methods and parameters within a given batch. Doing a `snakemake` dry-run (using the `-n` flag) will show all final files to be generated. A dictionary of parameters for a given method and id number will be generated in the results folder, named `param_dict.csv`. This should be saved  if needed, since it will be overwritten with each snakemake run.

## Running in the terminal
The pipeline can be run using the `snakemake` command. For example:

```snakemake --cores all --use-conda```

will run the pipeline using all provided cores and create the conda environments for each job. It is a good idea to add the `-n` flag for a dry run to get a preview of what will actually be run.

## Cluster
To use this pipeline on a cluster it is recommended to use a cluster profile. An example for a SLURM cluster is in `configs/cluster/config.yaml`. To use this profile with `sbatch` one could use the following in a `.sh` file: 

```
#!/bin/bash

#SBATCH -o logs/snakemake_output_%j.txt
#SBATCH -e logs/snakemake_output_%j.txt
#SBATCH -J tx_snakemake
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3-00:00:00
#SBATCH --nice=10000
#SBATCH -p cpu_p

mkdir logs

conda activate tx

snakemake --profile configs/cluster/
```

Here, the `snakemake` command will be run on a compute node so it does not terminate if the user closes the submit node. The snakemake job will then in turn submit more jobs to the cluster using the arguments specified in the cluster profile. In this case, the `snakemake` command itself will use the sbatch arguments, i.e. a time limit of 3 days. However, the sbatch arguments for each job will be specified by the `Snakefile` and the cluster `config.yaml` (i.e. the job `watershed` will have a time limit of 12 hours). 

Tip for cluster users: if the snakemake job unexpectedly terminates, make sure to unlock the working directory before re-running the command.

# Input and Intermediate files

The pipeline will generate many intermediate files in the `RESULTS` folder specified in the config. Files will follow this naming convention:
`[type]_[method1]-[id]_[method2]-[id].[ext]`

- `type`: the kind of data in the file, such as "segments" or "counts"
- `method2`: the method used to generate that file, such as `cellpose`. This will be the last method listed
- `method1`: these are the previous methods used to generate that data. For example, an `assignments.csv` file will probably have come from a `segments.tif` file and its associated methods
- `id`: this is the specific id code for a method that details its parameters. This information will be in the `param_dict.csv` made upon running the `Snakefile`

## Format of input files

- `image`: The DAPI stain of the tissue as a single channel `tif`
- `sc_data`: The scRNA-seq data as an AnnData `.h5ad` file with the cell type in the key `adata.obs['celltype']`
- `molecules`: The spatial transcriptomics data as a `.csv` with 4 columns in this order: 1) Gene name, 2) X-coord, 3) Y-coord, 4) Cell type. The x and y coordinates should be in pixel coordinates of the DAPI image.
- `segmented_image`: Used with the `custom` segmentation key- should be a segmented image with background = 0. If the cells are binary they should be seperated (unique values will be assigned during the `segmentation` step).

## Format of intermediate files

- `segments.tif`: a `.tif` representing segmented cells as a matrix. Each value in the matrix corresponds to a unique cell with 0 being background(Note: does not always start from 1). This is an overlay of the original DAPI stain.
- `areas.csv`: the area of each cell from the `segments.tif` in pixels
- `assignments.csv`: similar to `molecules.csv` but now with a column for the id of the cell each molecule is assigned to
- `celltypes.csv`: the type for each cell given an assignment algorithm- not generated by every assignment method
- `counts.h5ad`: The normalized AnnData count matrix for the spatial data. Raw counts stored in `adata.layers['raw_counts']`. Cell area stored in `adata.obs['area']`, cell type stored in `adata.obs['celltype']` and `adata.obs['prior_celltype']`
- `metrics.csv`: Selected metrics calculated for a given `counts.h5ad`
