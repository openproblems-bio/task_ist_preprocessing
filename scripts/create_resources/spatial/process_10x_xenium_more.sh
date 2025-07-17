#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "10x_xenium/2023_10x_human_lung_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE/Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_outs.zip
    dataset_name: "Xenium Preview Human Non diseased Lung With Add on FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human lung sections, using a development version of the Human Lung Gene Expression Panel and custom add-on panel for lung cell types."
    dataset_description: "Preview of Xenium In Situ Gene Expression data for adult human lung sections, using a development version of the Human Lung Gene Expression Panel and custom add-on panel for lung cell types. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.3.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_lung_cancer_xenium"
    input: https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.3.0/Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE/Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE_outs.zip
    dataset_name: "Xenium Preview Human Lung Cancer With Add on 2 FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human lung sections, using a development version of the Human Lung Gene Expression Panel and custom add-on panel for lung cell types."
    dataset_description: "Preview of Xenium In Situ Gene Expression data for adult human lung sections, using a development version of the Human Lung Gene Expression Panel and custom add-on panel for lung cell types. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.3.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_pancreas_cancer_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.6.0/Xenium_V1_hPancreas_Cancer_Add_on_FFPE/Xenium_V1_hPancreas_Cancer_Add_on_FFPE_outs.zip
    dataset_name: "Xenium V1 hPancreas Cancer Add on FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/pancreatic-cancer-with-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for human pancreatic cancer sections using the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Xenium In Situ Gene Expression data for human pancreatic cancer sections using the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.6.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_brain_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon_outs.zip
    dataset_name: "Xenium V1 FFPE Human Brain Healthy With Addon"
    dataset_url: "https://www.10xgenomics.com/datasets/xenium-human-brain-preview-data-1-standard"
    dataset_summary: "Adult human brain cortical section, healthy."
    dataset_description: "Preview of Xenium In Situ Gene Expression data for adult human brain cortical sections, using a development version of the Human Brain Gene Expression Panel and two custom add-on panels for brain cell types. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.3.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_human_skin_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hSkin_nondiseased_section_1_FFPE/Xenium_V1_hSkin_nondiseased_section_1_FFPE_outs.zip
    dataset_name: "Xenium V1 hSkin nondiseased section 1 FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-skin-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for adult human skin sections using the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Xenium In Situ Gene Expression data for adult human skin sections using the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.9.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_human_liver_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hLiver_nondiseased_section_FFPE/Xenium_V1_hLiver_nondiseased_section_FFPE_outs.zip
    dataset_name: "Xenium V1 hLiver nondiseased section FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for adult human liver sections using the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Xenium In Situ Gene Expression data for adult human liver sections using the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.9.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_human_liver_cancer_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hLiver_cancer_section_FFPE/Xenium_V1_hLiver_cancer_section_FFPE_outs.zip
    dataset_name: "Xenium V1 hLiver cancer section FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for adult human liver sections using the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Xenium In Situ Gene Expression data for adult human liver sections using the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.9.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_human_heart_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hHeart_nondiseased_section_FFPE/Xenium_V1_hHeart_nondiseased_section_FFPE_outs.zip
    dataset_name: "Xenium V1 hHeart nondiseased section FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-heart-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for adult human heart section using the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Xenium In Situ Gene Expression data for adult human heart section using the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.9.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_colon_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.6.0/Xenium_V1_hColon_Non_diseased_Add_on_FFPE/Xenium_V1_hColon_Non_diseased_Add_on_FFPE_outs.zip
    dataset_name: "Xenium V1 hColon Non diseased Add on FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-colon-preview-data-xenium-human-colon-gene-expression-panel-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human colon sections, using a development version of the Xenium Human Colon Gene Expression Panel."
    dataset_description: "Preview of Xenium In Situ Gene Expression data for adult human colon sections, using a development version of the Xenium Human Colon Gene Expression Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.6.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_colon_cancer_xenium"
    input: https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.6.0/Xenium_V1_hColon_Cancer_Add_on_FFPE/Xenium_V1_hColon_Cancer_Add_on_FFPE_outs.zip
    dataset_name: "Xenium V1 hColon Cancer Add on FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/human-colon-preview-data-xenium-human-colon-gene-expression-panel-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human colon sections, using a development version of the Xenium Human Colon Gene Expression Panel."
    dataset_description: "Cancer; stage 2A adenocarcinoma. Preview of Xenium In Situ Gene Expression data for adult human colon sections, using a development version of the Xenium Human Colon Gene Expression Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.6.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_kidney_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hKidney_nondiseased_section/Xenium_V1_hKidney_nondiseased_section_outs.zip
    dataset_name: "Xenium V1 hKidney nondiseased section"
    dataset_url: "https://www.10xgenomics.com/datasets/human-kidney-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human kidney sections, using a development version of the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Preview of Xenium In Situ Gene Expression data for adult human kidney sections, using a development version of the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.5.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_kidney_cancer_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hKidney_cancer_section/Xenium_V1_hKidney_cancer_section_outs.zip
    dataset_name: "Xenium V1 hKidney cancer section"
    dataset_url: "https://www.10xgenomics.com/datasets/human-kidney-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human kidney sections, using a development version of the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "Kidney cancer (papillary renal cell carcinoma, or PRCC). Preview of Xenium In Situ Gene Expression data for adult human kidney sections, using a development version of the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.5.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_mouse_colon_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_mouse_Colon_FF/Xenium_V1_mouse_Colon_FF_outs.zip
    dataset_name: "Xenium V1 mouse Colon FF"
    dataset_url: "https://www.10xgenomics.com/datasets/fresh-frozen-mouse-colon-with-xenium-multimodal-cell-segmentation-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression with Cell Segmentation Staining data for mouse colon tissue using the Xenium Mouse Tissue Atlassing Panel."
    dataset_description: "Xenium In Situ Gene Expression with Cell Segmentation Staining data for mouse colon tissue using the Xenium Mouse Tissue Atlassing Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 2.0.0"
    dataset_organism: "mus_musculus"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2023_10x_human_lymph_node_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hLymphNode_nondiseased_section/Xenium_V1_hLymphNode_nondiseased_section_outs.zip
    dataset_name: "Xenium V1 hLymphNode nondiseased section"
    dataset_url: "https://www.10xgenomics.com/datasets/human-lymph-node-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
    dataset_summary: "Preview of Xenium In Situ Gene Expression data for adult human lymph node, using a development version of the Xenium Human Multi-Tissue and Cancer Panel."
    dataset_description: "The selected section represents a non-diseased lymph node tissue. Preview of Xenium In Situ Gene Expression data for adult human lymph node, using a development version of the Xenium Human Multi-Tissue and Cancer Panel. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 1.5.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]

  - id: "10x_xenium/2024_10x_human_ovarian_cancer_xenium"
    input: https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Ovarian_Cancer_Addon_FFPE/Xenium_V1_Human_Ovarian_Cancer_Addon_FFPE_outs.zip
    dataset_name: "Xenium V1 Human Ovarian Cancer Addon FFPE"
    dataset_url: "https://www.10xgenomics.com/datasets/ffpe-human-ovarian-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard"
    dataset_summary: "Xenium In Situ Gene Expression data for human ovarian cancer tissue using the Xenium Human Immuno-Oncology Profiling Panel with custom add-on."
    dataset_description: "Ovary Serous Carcinoma; Stage II-A; Grade 3. Xenium In Situ Gene Expression data for human ovarian cancer tissue using the Xenium Human Immuno-Oncology Profiling Panel with custom add-on. In Situ Gene Expression dataset analyzed using Xenium Onboard Analysis 2.0.0"
    dataset_organism: "homo_sapiens"
    segmentation_id: [cell, nucleus]


output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_tenx_xenium/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,10x_xenium
