
To run the vizgen merscope workflow (`src/datasets/workflows/process_vizgen_merscope/test.sh`) a google cloud authentication is needed. See [here](https://www.nextflow.io/docs/latest/google.html#credentials) for more details. In short, run the command:
```
gcloud auth application-default login
``` 

Note that an install of gcloud-cli in a conda environment is not sufficient for nextflow to point to the correct authentication files. Instead a system wide installation is needed. E.g. for mac `brew install --cask gcloud-cli`.