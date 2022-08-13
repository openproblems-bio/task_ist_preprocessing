# Data Sources
## Developmental Heart Data:
[Asp et al. 2019, Cell](https://www.sciencedirect.com/science/article/pii/S0092867419312826) 

[Images](https://figshare.com/articles/dataset/ISS_data_in_A_spatiotemporal_organ-wide_gene_expression_and_cell_atlas_of_the_developing_human_heart_/10058048/1)

[Other Data](https://data.mendeley.com/datasets/mbvhhf8m62/2)

## Mouse Hippocampus:
[Qian et al. 2020, Nature Methods](https://www.nature.com/articles/s41592-019-0631-4)

[Images](https://figshare.com/articles/dataset/CA1_neuron_cell_typing_results_of_pciSeq_probabilistic_cell_typing_by_in_situ_sequencing_/7150760/1)

[Other Data](https://su.figshare.com/articles/dataset/pciSeq_files_in_csv_format/10318610/1)

### Notes
For hippocampus, images/spot were downloaded directly. The single cell data was opened in MATLAB using their [repo](https://github.com/cortex-lab/Transcriptomics.git). Each field of the `GeneSet` object was saved seperately and reloaded via a Jupyter notebook into an `AnnData` object. There were duplicate observation and variable names. The observations were made unique via `obs_names_make_unique()` from `AnnData`. None of the duplicate genes (`AnnData.var`) appeared in the spatial data, so they were also made unique in the same way.

## Fetal Liver/Hemapoietic Stem Cells:
[Lu et al. 2021, Cell Discovery](https://www.nature.com/articles/s41421-021-00266-1)

[Data](https://campuspress.yale.edu/wanglab/HSCMERFISH/)

## Misc
For testing purposes a small, cropped section of the developmental heart data (PCW4.5_1) was occasionally used 
