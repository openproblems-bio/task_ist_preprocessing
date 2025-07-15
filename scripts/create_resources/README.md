


The spatial datasets can be downloaded and processed with scripts in `task_ist_preprocessing/scripts/create_resources/spatial` 

The scRNAseq reference datasets can be downloaded and processed with scripts in `task_ist_preprocessing/scripts/create_resources/sc` 

Note that scripts in those two dirs create the datasets at:  
```
s3://openproblems-data/resources/datasets/*/
```

Finally the script `task_ist_preprocessing/scripts/create_resources/combine/process_datasets.sh` creates the combined, processed datasets (pairs of spatial and scRNAseq). This processing/combining step subsets the genes to the common genes and saves metadata of the combined object in the output files. The combined datasets are saved at 
```
s3://openproblems-data/resources/task_ist_preprocessing/datasets/
```
(The final benchmark runs on all datasets in this location)