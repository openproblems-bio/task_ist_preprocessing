# txsim-pipeline Development
Guide to helping develop the `txsim-pipeline`

# Adding New Method
This applies to a new method developed for segmentation, dot assignment, normalization etc. Examples are pciSeq and Baysor. 
The easiest types of methods to add are those that exist as python packages on [PyPI](https://pypi.org/) or as a [conda](https://docs.conda.io/en/latest/) package. Other methods either not written in python or those that are not packages, require a bit more effort to make them compatible with the pipeline.

## For all methods

### Snakemake
All methods will first require a new `snakemake` rule. Information on the general syntax of `snakemake` can be found [here](https://snakemake.readthedocs.io/en/stable/index.html). For this project each method has its own rule in the `Snakefile`, and should be named the same as the name of the method. The current convention for `txsim` is run a python script (see the `scripts` folder) using a shell command with command line arguments, and then parse said arguments in the python script (see below). The script is run in an environment specified in `envs`.

### Parameters
All methods should have their default parameters specified in the `defaults.yaml` in `configs`. These will not be passed into the method, and are for the readable `params` dictionary after a run. This is useful when sorting outputs.
Each method may have a unique way of passing in parameters, however, to work with the pipeline they must be formatted in `config.yaml` in `config`. The recommended way is to have the following hierarchy:

```
batch:
  dataset: [dataset_1]
  workflow: [group_1, group_2, group_3]
  group_1:
    method_1:
        parameter_1: 0
        parameter_2: True
        parameter_3: 2.5
    method_2:
        parameter_1: False
  group_2:
    method_3:
  group_3:
    method_4:
  group_3_params:
    parameter_1: True
    parameter_2: 5
 ...
```

With this structure, the `get_params` function within the `Snakefile` will return a method-specific parameter dictionary (using `'p'` as a third argument) or a group parameter for that method (using the parameter name as a third argument). For instance, `get_params` would return 
```
(dict) {parameter_1: 0, parameter_2: True, parameter_3: 2.5}
```
for `method_1` using the `'p'` argument. Likewise, `get_params` would return `5` using the `parameter_2` argument for `method_4` since those params apply to the whole group. This function is very useful for inputting parameters into a snakemake rule. 

### Scripts
Currently, each method has its own python script, or a few related methods will share one script. Within the script, specific methods are run, either using a wrapper in the `txsim` package,  an imported library, or through the command line interface. The last is used for non-python methods (i.e. Baysor). After the method is run, the outputs must be saved from the python script. Specifying the output in the `Snakefile` will not save the file automatically. The naming conventions for input and output are specified in the `README.md`. For a new method, it is recommended to add a new python script corresponding to that method.


## For python methods

### Adding a new package dependency 
Many methods are written as python packages on PyPI or conda. To utilize a new package/library, simply add a new conda environment in the `envs` folder, and use this environment in the `snakemake` rule. If adding a new package to integrate into the main `txsim` library, be sure to update the dependencies in `pyproject.toml` accordingly. 

### Adding a library not as a package. 
There are many ways to add a local library or repository, including locally importing the package via `pip` using the `-e` flag. The exact way of integrating a local library will depend on the library itself.

## For non-python methods
Adding a non-python method may be more difficult. Many `R` methods can be installed from `conda`, or even run in python. However, some methods may require the creation of a new container (i.e. a Docker container), which can be run through `snakemake` using `Singularity` (specified by the `container` flag in the rule). The exact implementation will vary greatly from method to method, but Baysor is an example of such a method. See the snakemake [documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#running-jobs-in-containers) for more info on containers.

# Editing other aspects of the pipeline
The pipeline itself mainly contains scripts and files specifying the workflow organization. In order to edit specific functions of `txsim` (such as how counts are normalized), directly edit the `txsim` package. Adding datasets is described in the main README. For the general structure of the repository:
- The `scripts` folder contains python scripts which run individual methods and parse command line arguments
- The `envs` folder contains environments used by `snakemake` to run the scripts
- The `configs` folder contains config files and defaults 
- `TxsimConfig.py` contains a ParsedConfig class which parses the config file, establishes the method combinations, and saves the parameters as readable csv files.
- `Snakefile` specifies the snakemake rules which are used to determine which scripts to run and which parameters to use.
