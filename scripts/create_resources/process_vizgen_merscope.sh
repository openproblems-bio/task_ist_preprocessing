#!/bin/bash

# Get the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

publish_dir="s3://openproblems-data/resources/datasets"

cat > /tmp/params.yaml << HERE
param_list:

  - id: "vizgen_merscope/2022_vizgen_human_breast_cancer_merfish/rep1"
    gcloud_bucket: "vz-ffpe-showcase"
    dataset_bucket_name: "HumanBreastCancerPatient1"
    dataset_name: "2022 Vizgen Human Breast Cancer MERFISH Patient1"
    dataset_url: "https://info.vizgen.com/ffpe-showcase?submissionGuid=a93dbab5-c128-4269-afe3-82ea2bf9cdaf"
    dataset_summary: "Human Breast Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
    dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. This data release includes 16 MERFISH datasets generated by the MERSCOPE Platform from 8 different human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "vizgen_merscope/2022_vizgen_human_liver_cancer_merfish/rep1"
    gcloud_bucket: "vz-ffpe-showcase"
    dataset_bucket_name: "HumanLiverCancerPatient1"
    dataset_name: "2022 Vizgen Human Liver Cancer MERFISH Patient1"
    dataset_url: "https://info.vizgen.com/ffpe-showcase?submissionGuid=a93dbab5-c128-4269-afe3-82ea2bf9cdaf"
    dataset_summary: "Human Liver Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
    dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. This data release includes 16 MERFISH datasets generated by the MERSCOPE Platform from 8 different human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "vizgen_merscope/2022_vizgen_human_liver_cancer_merfish/rep2"
    gcloud_bucket: "vz-ffpe-showcase"
    dataset_bucket_name: "HumanLiverCancerPatient2"
    dataset_name: "2022 Vizgen Human Liver Cancer MERFISH Patient2"
    dataset_url: "https://info.vizgen.com/ffpe-showcase?submissionGuid=a93dbab5-c128-4269-afe3-82ea2bf9cdaf"
    dataset_summary: "Human Liver Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
    dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. This data release includes 16 MERFISH datasets generated by the MERSCOPE Platform from 8 different human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "vizgen_merscope/2022_vizgen_human_lung_cancer_merfish/rep1"
    gcloud_bucket: "vz-ffpe-showcase"
    dataset_bucket_name: "HumanLungCancerPatient1"
    dataset_name: "2022 Vizgen Human Lung Cancer MERFISH Patient1"
    dataset_url: "https://info.vizgen.com/ffpe-showcase?submissionGuid=a93dbab5-c128-4269-afe3-82ea2bf9cdaf"
    dataset_summary: "Human Lung Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
    dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. This data release includes 16 MERFISH datasets generated by the MERSCOPE Platform from 8 different human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

  - id: "vizgen_merscope/2022_vizgen_human_lung_cancer_merfish/rep2"
    gcloud_bucket: "vz-ffpe-showcase"
    dataset_bucket_name: "HumanLungCancerPatient2"
    dataset_name: "2022 Vizgen Human Lung Cancer MERFISH Patient2"
    dataset_url: "https://info.vizgen.com/ffpe-showcase?submissionGuid=a93dbab5-c128-4269-afe3-82ea2bf9cdaf"
    dataset_summary: "Human Lung Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
    dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. This data release includes 16 MERFISH datasets generated by the MERSCOPE Platform from 8 different human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
    dataset_organism: "homo_sapiens"
    segmentation_id: ["cell"]

output_dataset: "\$id/dataset.zarr"
output_state: "\$id/state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_ist_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_vizgen_merscope/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels datasets,vizgen_merscope



# More datasets that can be simply added:
# TODO: Make a decision on replicate naming (see ovarian cancer replicate that has multiple slices)

#   - id: "vizgen_merscope/2022_vizgen_human_colon_cancer_merfish/rep1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanColonCancerPatient1"
#     dataset_name: "2022 Vizgen Human Colon Cancer MERFISH Patient1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Colon Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "The MERSCOPE FFPE Human Immuno-Oncology Data Release was generated using the MERSCOPE FFPE Sample Prep Solution and the MERSCOPE Immuno-Oncology Predesigned Panel. It includes datasets from various human tumor types, each measuring 500 genes representing approximately 4 billion transcripts and 9 million cells cumulatively."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_colon_cancer_merfish/rep2"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanColonCancerPatient2"
#     dataset_name: "2022 Vizgen Human Colon Cancer MERFISH Patient2"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Colon Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_melanoma_merfish/rep1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanMelanomaPatient1"
#     dataset_name: "2022 Vizgen Human Melanoma MERFISH Patient1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Melanoma data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_melanoma_merfish/rep2"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanMelanomaPatient2"
#     dataset_name: "2022 Vizgen Human Melanoma MERFISH Patient2"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Melanoma data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_ovarian_cancer_merfish/rep1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanOvarianCancerPatient1"
#     dataset_name: "2022 Vizgen Human Ovarian Cancer MERFISH Patient1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Ovarian Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   # Patient 2 has multiple slices
#   - id: "vizgen_merscope/2022_vizgen_human_ovarian_cancer_merfish/rep2_slice1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanOvarianCancerPatient2Slice1"
#     dataset_name: "2022 Vizgen Human Ovarian Cancer MERFISH Patient2 Slice1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Ovarian Cancer data (Slice 1) from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_ovarian_cancer_merfish/rep2_slice2"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanOvarianCancerPatient2Slice2"
#     dataset_name: "2022 Vizgen Human Ovarian Cancer MERFISH Patient2 Slice2"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Ovarian Cancer data (Slice 2) from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_ovarian_cancer_merfish/rep2_slice3"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanOvarianCancerPatient2Slice3"
#     dataset_name: "2022 Vizgen Human Ovarian Cancer MERFISH Patient2 Slice3"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Ovarian Cancer data (Slice 3) from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_prostate_cancer_merfish/rep1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanProstateCancerPatient1"
#     dataset_name: "2022 Vizgen Human Prostate Cancer MERFISH Patient1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Prostate Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_prostate_cancer_merfish/rep2"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanProstateCancerPatient2"
#     dataset_name: "2022 Vizgen Human Prostate Cancer MERFISH Patient2"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Prostate Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_uterine_cancer_merfish/rep1"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanUterineCancerPatient1"
#     dataset_name: "2022 Vizgen Human Uterine Cancer MERFISH Patient1"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Uterine Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]

#   - id: "vizgen_merscope/2022_vizgen_human_uterine_cancer_merfish/rep2"
#     gcloud_bucket: "vz-ffpe-showcase"
#     dataset_bucket_name: "HumanUterineCancerPatient2"
#     dataset_name: "2022 Vizgen Human Uterine Cancer MERFISH Patient2"
#     dataset_url: "https://info.vizgen.com/ffpe-showcase"
#     dataset_summary: "Human Uterine Cancer data from the MERSCOPE FFPE Human Immuno-Oncology Data Release."
#     dataset_description: "Same as above."
#     dataset_organism: "homo_sapiens"
#     segmentation_id: ["cell"]
