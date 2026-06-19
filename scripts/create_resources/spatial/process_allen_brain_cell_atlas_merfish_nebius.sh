#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "allen_brain_cell_atlas_merfish/mouse1_coronal/rep1"
    mouse: "mouse1_coronal"
    experiment_id: "220422_wb3_co1_1B_6z18R_merfish5"
    dataset_name: "Allen Brain Cell Atlas MERFISH Mouse 1 Coronal"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_summary: "Brain-wide MERFISH spatial transcriptomics of mouse 1 (coronal), Allen Brain Cell Atlas."
    dataset_description: "Brain-wide MERFISH spatial transcriptomics data from the Zhuang lab. Mouse 1 coronal section imaged with ~1100 gene panel."
    dataset_organism: "mus_musculus"
    dataset_reference: "@article{Yao2023, author={Yao, Zizhen and others}, title={A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain}, journal={Nature}, year={2023}}"
    segmentation_id: ["cell"]

  - id: "allen_brain_cell_atlas_merfish/mouse2_coronal/rep1"
    mouse: "mouse2_coronal"
    experiment_id: "220601_wb3_co2_1_5z18R2bd_merfish5"
    dataset_name: "Allen Brain Cell Atlas MERFISH Mouse 2 Coronal"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_summary: "Brain-wide MERFISH spatial transcriptomics of mouse 2 (coronal), Allen Brain Cell Atlas."
    dataset_description: "Brain-wide MERFISH spatial transcriptomics data from the Zhuang lab. Mouse 2 coronal section imaged with ~1100 gene panel."
    dataset_organism: "mus_musculus"
    dataset_reference: "@article{Yao2023, author={Yao, Zizhen and others}, title={A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain}, journal={Nature}, year={2023}}"
    segmentation_id: ["cell"]

  - id: "allen_brain_cell_atlas_merfish/mouse3_sagittal/rep1"
    mouse: "mouse3_sagittal"
    experiment_id: "220609_wb3_sa1_1_5z18R_merfish5"
    dataset_name: "Allen Brain Cell Atlas MERFISH Mouse 3 Sagittal"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_summary: "Brain-wide MERFISH spatial transcriptomics of mouse 3 (sagittal), Allen Brain Cell Atlas."
    dataset_description: "Brain-wide MERFISH spatial transcriptomics data from the Zhuang lab. Mouse 3 sagittal section imaged with ~1100 gene panel."
    dataset_organism: "mus_musculus"
    dataset_reference: "@article{Yao2023, author={Yao, Zizhen and others}, title={A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain}, journal={Nature}, year={2023}}"
    segmentation_id: ["cell"]

  - id: "allen_brain_cell_atlas_merfish/mouse4_sagittal/rep1"
    mouse: "mouse4_sagittal"
    experiment_id: "220912_wb3_sa2_2_5z18R_merfish5"
    dataset_name: "Allen Brain Cell Atlas MERFISH Mouse 4 Sagittal"
    dataset_url: "https://download.brainimagelibrary.org/29/3c/293cc39ceea87f6d/"
    dataset_summary: "Brain-wide MERFISH spatial transcriptomics of mouse 4 (sagittal), Allen Brain Cell Atlas."
    dataset_description: "Brain-wide MERFISH spatial transcriptomics data from the Zhuang lab. Mouse 4 sagittal section imaged with ~1100 gene panel."
    dataset_organism: "mus_musculus"
    dataset_reference: "@article{Yao2023, author={Yao, Zizhen and others}, title={A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain}, journal={Nature}, year={2023}}"
    segmentation_id: ["cell"]

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_allen_brain_cell_atlas_merfish/main.nf \
  --workspace 167877437119966 \
  --compute-env 5hfmdCBxMRd4nHZaJKYEQZ \
  --params-file /tmp/params.yaml \
  --config src/base/labels_nebius.config \
  --labels datasets,allen_brain_cell_atlas_merfish
