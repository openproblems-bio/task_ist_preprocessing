#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

input_dir="s3://openproblems-data/resources/datasets"
publish_dir="s3://openproblems-data/resources/task_ist_preprocessing/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "mouse_brain_combined/rep1"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep1/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_name: "Mouse brain combined 2023 tenx Xenium replicate 1 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;10.1038/s41586-023-06812-z"
    dataset_summary: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1);A high-resolution scRNAseq atlas of cell types in the whole mouse brain"
    dataset_description: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1). Replicate results demonstrate the high reproducibility of data generated by the platform. 10x Genomics obtained tissue from a C57BL/6 mouse from Charles River Laboratories. Three adjacent 10µm sections were placed on the same slide. Tissues were prepared following the demonstrated protocols Xenium In Situ for Fresh Frozen Tissues - Tissue Preparation Guide (CG000579) and Xenium In Situ for Fresh Frozen Tissues - Fixation & Permeabilization (CG000581).;See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset."
    dataset_organism: "mus_musculus"
  - id: "mouse_brain_combined/rep2"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep2/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_name: "Mouse brain combined 2023 tenx Xenium replicate 2 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;10.1038/s41586-023-06812-z"
    dataset_summary: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1);A high-resolution scRNAseq atlas of cell types in the whole mouse brain"
    dataset_description: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1). Replicate results demonstrate the high reproducibility of data generated by the platform. 10x Genomics obtained tissue from a C57BL/6 mouse from Charles River Laboratories. Three adjacent 10µm sections were placed on the same slide. Tissues were prepared following the demonstrated protocols Xenium In Situ for Fresh Frozen Tissues - Tissue Preparation Guide (CG000579) and Xenium In Situ for Fresh Frozen Tissues - Fixation & Permeabilization (CG000581).;See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset."
    dataset_organism: "mus_musculus"
  - id: "mouse_brain_combined/rep3"
    input_sp: "$input_dir/10x_xenium/2023_10x_mouse_brain_xenium/rep3/dataset.zarr"
    input_sc: "$input_dir/allen_brain_cell_atlas/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"
    dataset_name: "Mouse brain combined 2023 tenx Xenium replicate 3 2023 Yao scRNAseq"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717"
    dataset_reference: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard;10.1038/s41586-023-06812-z"
    dataset_summary: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1);A high-resolution scRNAseq atlas of cell types in the whole mouse brain"
    dataset_description: "Demonstration of gene expression profiling for fresh frozen mouse brain on the Xenium platform using the pre-designed Mouse Brain Gene Expression Panel (v1). Replicate results demonstrate the high reproducibility of data generated by the platform. 10x Genomics obtained tissue from a C57BL/6 mouse from Charles River Laboratories. Three adjacent 10µm sections were placed on the same slide. Tissues were prepared following the demonstrated protocols Xenium In Situ for Fresh Frozen Tissues - Tissue Preparation Guide (CG000579) and Xenium In Situ for Fresh Frozen Tissues - Fixation & Permeabilization (CG000581).;See dataset_reference for more information. Note that we only took the 10xv2 data from the dataset."
    dataset_organism: "mus_musculus"

output_sc: "\$id/output_sc.h5ad"
output_sp: "\$id/output_sp.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_ist_preprocessing,process_datasets
