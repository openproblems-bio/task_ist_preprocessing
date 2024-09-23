# iST Preprocessing


<!--
This file is automatically generated from the tasks's api/*.yaml files.
Do not edit this file directly.
-->

Benchmarking approaches for preprocessing imaging-based spatial
transcriptomics

Repository:
[openproblems-bio/task_ist_preprocessing](https://github.com/openproblems-bio/task_ist_preprocessing)

## Description

Provide a clear and concise description of your task, detailing the
specific problem it aims to solve. Outline the input data types, the
expected output, and any assumptions or constraints. Be sure to explain
any terminology or concepts that are essential for understanding the
task.

Explain the motivation behind your proposed task. Describe the
biological or computational problem you aim to address and why it’s
important. Discuss the current state of research in this area and any
gaps or challenges that your task could help address. This section
should convince readers of the significance and relevance of your task.

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Louis Kümmerle    | author, maintainer |
| Malte D. Luecken  | author             |
| Daniel Strobl     | author             |
| Robrecht Cannoodt | author             |

## API

``` mermaid
flowchart TB
  file_common_ist("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-common-ist-dataset'>Common iST Dataset</a>")
  comp_data_preprocessor[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-data-preprocessor'>Data preprocessor</a>"/]
  file_raw_ist("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-raw-ist-dataset'>Raw iST Dataset</a>")
  file_scrnaseq_reference("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-scrna-seq-reference'>scRNA-seq Reference</a>")
  comp_method_segmentation[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-segmentation'>Segmentation</a>"/]
  comp_method_transcript_assignment[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-assignment'>Assignment</a>"/]
  comp_method_cell_type_annotation[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-cell-type-annotation'>Cell Type Annotation</a>"/]
  comp_method_expression_correction[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-expression-correction'>Expression correction</a>"/]
  comp_metric_similarity[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-metric'>Metric</a>"/]
  file_segmentation("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-segmentation'>Segmentation</a>")
  file_transcript_assignments("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-transcript-assignment'>Transcript Assignment</a>")
  file_spatial_with_cell_types("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-spatial-with-cell-types'>Spatial with Cell Types</a>")
  file_spatial_corrected_counts("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-spatial-corrected'>Spatial Corrected</a>")
  file_score("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-score'>Score</a>")
  comp_method_calculate_cell_volume[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-calculate-cell-volume'>Calculate Cell Volume</a>"/]
  comp_method_count_aggregation[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-count-aggregation'>Count Aggregation</a>"/]
  comp_metric_quality[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-quality-metric'>Quality Metric</a>"/]
  file_cell_volumes("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-cell-volumes'>Cell Volumes</a>")
  file_spatial_aggregated_counts("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-aggregated-counts'>Aggregated Counts</a>")
  comp_method_normalization[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-normalization'>Normalization</a>"/]
  comp_method_qc_filter[/"<a href='https://github.com/openproblems-bio/task_ist_preprocessing#component-type-qc-filter'>QC Filter</a>"/]
  file_spatial_normalized_counts("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-spatial-normalized'>Spatial Normalized</a>")
  file_spatial_qc_col("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-qc-columns'>QC Columns</a>")
  file_common_scrnaseq("<a href='https://github.com/openproblems-bio/task_ist_preprocessing#file-format-common-sc-dataset'>Common SC Dataset</a>")
  file_common_ist---comp_data_preprocessor
  comp_data_preprocessor-->file_raw_ist
  comp_data_preprocessor-->file_scrnaseq_reference
  file_raw_ist---comp_method_segmentation
  file_raw_ist---comp_method_transcript_assignment
  file_scrnaseq_reference-.-comp_method_transcript_assignment
  file_scrnaseq_reference-.-comp_method_cell_type_annotation
  file_scrnaseq_reference-.-comp_method_expression_correction
  file_scrnaseq_reference---comp_metric_similarity
  comp_method_segmentation-->file_segmentation
  comp_method_transcript_assignment-->file_transcript_assignments
  comp_method_cell_type_annotation-->file_spatial_with_cell_types
  comp_method_expression_correction-->file_spatial_corrected_counts
  comp_metric_similarity-->file_score
  file_segmentation-.-comp_method_transcript_assignment
  file_transcript_assignments-.-comp_method_cell_type_annotation
  file_transcript_assignments---comp_method_calculate_cell_volume
  file_transcript_assignments---comp_method_count_aggregation
  file_transcript_assignments---comp_metric_quality
  file_spatial_with_cell_types---comp_method_expression_correction
  file_spatial_corrected_counts---comp_metric_similarity
  file_spatial_corrected_counts---comp_metric_quality
  comp_method_calculate_cell_volume-->file_cell_volumes
  comp_method_count_aggregation-->file_spatial_aggregated_counts
  comp_metric_quality-->file_score
  file_cell_volumes-.-comp_method_normalization
  file_spatial_aggregated_counts---comp_method_normalization
  file_spatial_aggregated_counts---comp_method_qc_filter
  comp_method_normalization-->file_spatial_normalized_counts
  comp_method_qc_filter-->file_spatial_qc_col
  file_spatial_normalized_counts---comp_method_cell_type_annotation
  file_spatial_qc_col---comp_metric_similarity
  file_spatial_qc_col---comp_metric_quality
  file_common_scrnaseq---comp_data_preprocessor
```

## File format: Common iST Dataset

An unprocessed spatial imaging dataset stored as a zarr file.

Example file:
`resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr`

Description:

This dataset contains raw images, labels, points, shapes, and tables as
output by a dataset loader.

Format:

<div class="small">

</div>

Data structure:

<div class="small">

</div>

## Component type: Data preprocessor

Preprocess a common dataset for the benchmark.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_ist` | `file` | An unprocessed spatial imaging dataset stored as a zarr file. |
| `--input_scrnaseq` | `file` | An unprocessed dataset as output by a dataset loader. |
| `--output_ist` | `file` | (*Output*) A spatial transcriptomics dataset, preprocessed for this benchmark. |
| `--output_scrnaseq` | `file` | (*Output*) A single-cell reference dataset, preprocessed for this benchmark. |

</div>

## File format: Raw iST Dataset

A spatial transcriptomics dataset, preprocessed for this benchmark.

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/raw_ist.zarr`

Description:

This dataset contains preprocessed images, labels, points, shapes, and
tables for spatial transcriptomics data.

Format:

<div class="small">

</div>

Data structure:

<div class="small">

</div>

## File format: scRNA-seq Reference

A single-cell reference dataset, preprocessed for this benchmark.

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad`

Description:

This dataset contains preprocessed counts and metadata for single-cell
RNA-seq data.

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'cell_type_level2', 'cell_type_level3', 'cell_type_level4', 'dataset_id', 'assay', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage', 'development_stage_ontology_term_id', 'disease', 'disease_ontology_term_id', 'donor_id', 'is_primary_data', 'organism', 'organism_ontology_term_id', 'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id', 'sex', 'sex_ontology_term_id', 'suspension_type', 'tissue', 'tissue_ontology_term_id', 'tissue_general', 'tissue_general_ontology_term_id', 'batch', 'soma_joinid'
     var: 'feature_id', 'feature_name', 'soma_joinid', 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_type"]` | `string` | Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level2"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level3"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level4"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["dataset_id"]` | `string` | (*Optional*) Identifier for the dataset from which the cell data is derived, useful for tracking and referencing purposes. |
| `obs["assay"]` | `string` | (*Optional*) Type of assay used to generate the cell data, indicating the methodology or technique employed. |
| `obs["assay_ontology_term_id"]` | `string` | (*Optional*) Experimental Factor Ontology (`EFO:`) term identifier for the assay, providing a standardized reference to the assay type. |
| `obs["cell_type_ontology_term_id"]` | `string` | (*Optional*) Cell Ontology (`CL:`) term identifier for the cell type, offering a standardized reference to the specific cell classification. |
| `obs["development_stage"]` | `string` | (*Optional*) Stage of development of the organism or tissue from which the cell is derived, indicating its maturity or developmental phase. |
| `obs["development_stage_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the developmental stage, providing a standardized reference to the organism’s developmental phase. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Developmental Stages (`HsapDv:`) ontology is used. If the organism is mouse (`organism_ontology_term_id == 'NCBITaxon:10090'`), then the Mouse Developmental Stages (`MmusDv:`) ontology is used. Otherwise, the Uberon (`UBERON:`) ontology is used. |
| `obs["disease"]` | `string` | (*Optional*) Information on any disease or pathological condition associated with the cell or donor. |
| `obs["disease_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the disease, enabling standardized disease classification and referencing. Must be a term from the Mondo Disease Ontology (`MONDO:`) ontology term, or `PATO:0000461` from the Phenotype And Trait Ontology (`PATO:`). |
| `obs["donor_id"]` | `string` | (*Optional*) Identifier for the donor from whom the cell sample is obtained. |
| `obs["is_primary_data"]` | `boolean` | (*Optional*) Indicates whether the data is primary (directly obtained from experiments) or has been computationally derived from other primary data. |
| `obs["organism"]` | `string` | (*Optional*) Organism from which the cell sample is obtained. |
| `obs["organism_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the organism, providing a standardized reference for the organism. Must be a term from the NCBI Taxonomy Ontology (`NCBITaxon:`) which is a child of `NCBITaxon:33208`. |
| `obs["self_reported_ethnicity"]` | `string` | (*Optional*) Ethnicity of the donor as self-reported, relevant for studies considering genetic diversity and population-specific traits. |
| `obs["self_reported_ethnicity_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the self-reported ethnicity, providing a standardized reference for ethnic classifications. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Ancestry Ontology (`HANCESTRO:`) is used. |
| `obs["sex"]` | `string` | (*Optional*) Biological sex of the donor or source organism, crucial for studies involving sex-specific traits or conditions. |
| `obs["sex_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the biological sex, ensuring standardized classification of sex. Only `PATO:0000383`, `PATO:0000384` and `PATO:0001340` are allowed. |
| `obs["suspension_type"]` | `string` | (*Optional*) Type of suspension or medium in which the cells were stored or processed, important for understanding cell handling and conditions. |
| `obs["tissue"]` | `string` | (*Optional*) Specific tissue from which the cells were derived, key for context and specificity in cell studies. |
| `obs["tissue_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the tissue, providing a standardized reference for the tissue type. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`. |
| `obs["tissue_general"]` | `string` | (*Optional*) General category or classification of the tissue, useful for broader grouping and comparison of cell data. |
| `obs["tissue_general_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the general tissue category, aiding in standardizing and grouping tissue types. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`. |
| `obs["batch"]` | `string` | (*Optional*) A batch identifier. This label is very context-dependent and may be a combination of the tissue, assay, donor, etc. |
| `obs["soma_joinid"]` | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the cell. |
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature, usually a ENSEMBL gene id. |
| `var["feature_name"]` | `string` | A human-readable name for the feature, usually a gene symbol. |
| `var["soma_joinid"]` | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the feature. |
| `var["hvg"]` | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `var["hvg_score"]` | `double` | A score for the feature indicating how highly variable it is. |
| `obsm["X_pca"]` | `double` | The resulting PCA embedding. |
| `obsp["knn_distances"]` | `double` | K nearest neighbors distance matrix. |
| `obsp["knn_connectivities"]` | `double` | K nearest neighbors connectivities matrix. |
| `varm["pca_loadings"]` | `double` | The PCA loadings matrix. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["normalized"]` | `integer` | Normalized expression values. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. This is different from the `obs.dataset_id` field, which is the identifier for the dataset from which the cell data is derived. |
| `uns["dataset_name"]` | `string` | A human-readable name for the dataset. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## Component type: Segmentation

A segmentation of the spatial data into cells

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | A spatial transcriptomics dataset, preprocessed for this benchmark. |
| `--output` | `file` | (*Output*) A segmentation of a spatial transcriptomics dataset. |

</div>

## Component type: Assignment

Assigning transcripts to cells

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_ist` | `file` | A spatial transcriptomics dataset, preprocessed for this benchmark. |
| `--input_segmentation` | `file` | (*Optional*) A segmentation of a spatial transcriptomics dataset. |
| `--input_scrnaseq` | `file` | (*Optional*) A single-cell reference dataset, preprocessed for this benchmark. |
| `--output` | `file` | (*Output*) A spatial transcriptomics dataset with assigned transcripts. |

</div>

## Component type: Cell Type Annotation

Annotating cell types in spatial data

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_spatial_normalized_counts` | `file` | Normalized counts. |
| `--input_transcript_assignments` | `file` | (*Optional*) A spatial transcriptomics dataset with assigned transcripts. |
| `--input_scrnaseq_reference` | `file` | (*Optional*) A single-cell reference dataset, preprocessed for this benchmark. |
| `--celltype_key` | `string` | (*Optional*) NA. Default: `cell_type`. |
| `--output` | `file` | (*Output*) Normalized counts with cell type annotations. |

</div>

## Component type: Expression correction

Correcting expression levels in spatial data

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_spatial_with_cell_types` | `file` | Normalized counts with cell type annotations. |
| `--input_scrnaseq_reference` | `file` | (*Optional*) A single-cell reference dataset, preprocessed for this benchmark. |
| `--output` | `file` | (*Output*) Corrected spatial data counts with cell type annotations. |

</div>

## Component type: Metric

A metric for evaluating iST preprocessing methods

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | Corrected spatial data counts with cell type annotations. |
| `--input_qc_col` | `file` | QC columns for spatial data. |
| `--input_sc` | `file` | A single-cell reference dataset, preprocessed for this benchmark. |
| `--output` | `file` | (*Output*) Metric score file. |

</div>

## File format: Segmentation

A segmentation of a spatial transcriptomics dataset

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/segmentation.zarr`

Description:

This dataset contains a segmentation of the spatial transcriptomics
data.

Format:

<div class="small">

</div>

Data structure:

<div class="small">

</div>

## File format: Transcript Assignment

A spatial transcriptomics dataset with assigned transcripts

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/transcript_assignments.zarr`

Description:

This dataset contains the spatial transcriptomics data with assigned
transcripts.

Format:

<div class="small">

</div>

Data structure:

<div class="small">

</div>

## File format: Spatial with Cell Types

Normalized counts with cell type annotations

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_with_cell_types.h5ad`

Description:

This file contains the normalized counts of the spatial transcriptomics
data and cell type annotations.

Format:

<div class="small">

    AnnData object
     obs: 'cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'n_counts', 'n_genes', 'volume', 'cell_type'
     var: 'gene_name', 'n_counts', 'n_cells'
     layers: 'counts', 'normalized'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_id"]` | `string` | Unique identifier for the cell (from assignment step). |
| `obs["centroid_x"]` | `string` | X coordinate of the cell. |
| `obs["centroid_y"]` | `string` | Y coordinate of the cell. |
| `obs["centroid_z"]` | `string` | (*Optional*) Z coordinate of the cell. |
| `obs["n_counts"]` | `string` | Number of counts in the cell. |
| `obs["n_genes"]` | `string` | Number of genes in the cell. |
| `obs["volume"]` | `string` | Volume of the cell. |
| `obs["cell_type"]` | `string` | Cell type of the cell. |
| `var["gene_name"]` | `string` | Name of the gene. |
| `var["n_counts"]` | `string` | Number of counts of the gene. |
| `var["n_cells"]` | `string` | Number of cells expressing the gene. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["normalized"]` | `integer` | Normalized counts. |

</div>

## File format: Spatial Corrected

Corrected spatial data counts with cell type annotations

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_corrected_counts.h5ad`

Description:

This file contains the corrected counts of the spatial transcriptomics
data and cell type annotations.

Format:

<div class="small">

    AnnData object
     obs: 'cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'n_counts', 'n_genes', 'volume', 'cell_type'
     var: 'gene_name', 'n_counts', 'n_cells'
     layers: 'counts', 'normalized', 'normalized', 'normalized_uncorrected'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_id"]` | `string` | Unique identifier for the cell (from assignment step). |
| `obs["centroid_x"]` | `string` | X coordinate of the cell. |
| `obs["centroid_y"]` | `string` | Y coordinate of the cell. |
| `obs["centroid_z"]` | `string` | (*Optional*) Z coordinate of the cell. |
| `obs["n_counts"]` | `string` | Number of counts in the cell. |
| `obs["n_genes"]` | `string` | Number of genes in the cell. |
| `obs["volume"]` | `string` | Volume of the cell. |
| `obs["cell_type"]` | `string` | Cell type of the cell. |
| `var["gene_name"]` | `string` | Name of the gene. |
| `var["n_counts"]` | `string` | Number of counts of the gene. |
| `var["n_cells"]` | `string` | Number of cells expressing the gene. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["normalized"]` | `integer` | Normalized counts. |
| `layers["normalized"]` | `double` | (*Optional*) Corrected normalized expression. |
| `layers["normalized_uncorrected"]` | `double` | (*Optional*) Uncorrected normalized expression. |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'metric_ids', 'metric_values'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `uns["metric_ids"]` | `string` | One or more unique metric identifiers. |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

## Component type: Calculate Cell Volume

Calculate the volume of cells

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | A spatial transcriptomics dataset with assigned transcripts. |
| `--output` | `file` | (*Output*) An obs column of cell volumes calculated from spatial transcriptomics data. |

</div>

## Component type: Count Aggregation

Aggregating counts of transcripts within cells

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | A spatial transcriptomics dataset with assigned transcripts. |
| `--output` | `file` | (*Output*) Unprocessed raw counts after aggregation of transcripts to cells. |

</div>

## Component type: Quality Metric

A metric for evaluating the quality of the processed iST data

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | Corrected spatial data counts with cell type annotations. |
| `--input_qc_col` | `file` | QC columns for spatial data. |
| `--input_transcript_assignments` | `file` | A spatial transcriptomics dataset with assigned transcripts. |
| `--score` | `file` | (*Output*) Metric score file. |

</div>

## File format: Cell Volumes

An obs column of cell volumes calculated from spatial transcriptomics
data.

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/cell_volumes.h5ad`

Description:

An obs column of cell volumes calculated from spatial transcriptomics
data.

Format:

<div class="small">

    AnnData object
     obs: 'volume'

</div>

Data structure:

<div class="small">

| Slot            | Type     | Description             |
|:----------------|:---------|:------------------------|
| `obs["volume"]` | `string` | The volume of the cell. |

</div>

## File format: Aggregated Counts

Unprocessed raw counts after aggregation of transcripts to cells

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad`

Description:

This file contains the raw counts after aggregating transcripts to
cells.

Format:

<div class="small">

    AnnData object
     obs: 'cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'n_counts', 'n_genes'
     var: 'gene_name', 'n_counts', 'n_cells'
     layers: 'counts'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_id"]` | `string` | Unique identifier for the cell (from assignment step). |
| `obs["centroid_x"]` | `string` | X coordinate of the cell. |
| `obs["centroid_y"]` | `string` | Y coordinate of the cell. |
| `obs["centroid_z"]` | `string` | (*Optional*) Z coordinate of the cell. |
| `obs["n_counts"]` | `string` | Number of counts in the cell. |
| `obs["n_genes"]` | `string` | Number of genes in the cell. |
| `var["gene_name"]` | `string` | Name of the gene. |
| `var["n_counts"]` | `string` | Number of counts of the gene. |
| `var["n_cells"]` | `string` | Number of cells expressing the gene. |
| `layers["counts"]` | `integer` | Raw aggregated counts. |

</div>

## Component type: Normalization

Normalizing spatial transcriptomics data

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_spatial_aggregated_counts` | `file` | Unprocessed raw counts after aggregation of transcripts to cells. |
| `--input_cell_volumes` | `file` | (*Optional*) An obs column of cell volumes calculated from spatial transcriptomics data. |
| `--output` | `file` | (*Output*) Normalized counts. |

</div>

## Component type: QC Filter

Filtering cells based on QC metrics

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | Unprocessed raw counts after aggregation of transcripts to cells. |
| `--output` | `file` | (*Output*) QC columns for spatial data. |

</div>

## File format: Spatial Normalized

Normalized counts

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_normalized_counts.h5ad`

Description:

This file contains the normalized counts of the spatial transcriptomics
data.

Format:

<div class="small">

    AnnData object
     obs: 'cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'n_counts', 'n_genes'
     var: 'gene_name', 'n_counts', 'n_cells'
     layers: 'counts', 'normalized'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_id"]` | `string` | Unique identifier for the cell (from assignment step). |
| `obs["centroid_x"]` | `string` | X coordinate of the cell. |
| `obs["centroid_y"]` | `string` | Y coordinate of the cell. |
| `obs["centroid_z"]` | `string` | (*Optional*) Z coordinate of the cell. |
| `obs["n_counts"]` | `string` | Number of counts in the cell. |
| `obs["n_genes"]` | `string` | Number of genes in the cell. |
| `var["gene_name"]` | `string` | Name of the gene. |
| `var["n_counts"]` | `string` | Number of counts of the gene. |
| `var["n_cells"]` | `string` | Number of cells expressing the gene. |
| `layers["counts"]` | `integer` | Raw aggregated counts. |
| `layers["normalized"]` | `integer` | Normalized expression values. |

</div>

## File format: QC Columns

QC columns for spatial data

Example file:
`resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_qc_col.h5ad`

Description:

This file contains the QC-filter column for spatial data.

Format:

<div class="small">

    AnnData object
     obs: 'passed_QC'

</div>

Data structure:

<div class="small">

| Slot               | Type     | Description                                  |
|:-------------------|:---------|:---------------------------------------------|
| `obs["passed_QC"]` | `string` | Whether the cell passed the quality control. |

</div>

## File format: Common SC Dataset

An unprocessed dataset as output by a dataset loader.

Example file:
`resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad`

Description:

This dataset contains raw counts and metadata as output by a dataset
loader.

The format of this file is mainly derived from the [CELLxGENE schema
v4.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md).

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'cell_type_level2', 'cell_type_level3', 'cell_type_level4', 'dataset_id', 'assay', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage', 'development_stage_ontology_term_id', 'disease', 'disease_ontology_term_id', 'donor_id', 'is_primary_data', 'organism', 'organism_ontology_term_id', 'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id', 'sex', 'sex_ontology_term_id', 'suspension_type', 'tissue', 'tissue_ontology_term_id', 'tissue_general', 'tissue_general_ontology_term_id', 'batch', 'soma_joinid'
     var: 'feature_id', 'feature_name', 'soma_joinid', 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_type"]` | `string` | Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level2"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level3"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["cell_type_level4"]` | `string` | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `obs["dataset_id"]` | `string` | (*Optional*) Identifier for the dataset from which the cell data is derived, useful for tracking and referencing purposes. |
| `obs["assay"]` | `string` | (*Optional*) Type of assay used to generate the cell data, indicating the methodology or technique employed. |
| `obs["assay_ontology_term_id"]` | `string` | (*Optional*) Experimental Factor Ontology (`EFO:`) term identifier for the assay, providing a standardized reference to the assay type. |
| `obs["cell_type_ontology_term_id"]` | `string` | (*Optional*) Cell Ontology (`CL:`) term identifier for the cell type, offering a standardized reference to the specific cell classification. |
| `obs["development_stage"]` | `string` | (*Optional*) Stage of development of the organism or tissue from which the cell is derived, indicating its maturity or developmental phase. |
| `obs["development_stage_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the developmental stage, providing a standardized reference to the organism’s developmental phase. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Developmental Stages (`HsapDv:`) ontology is used. If the organism is mouse (`organism_ontology_term_id == 'NCBITaxon:10090'`), then the Mouse Developmental Stages (`MmusDv:`) ontology is used. Otherwise, the Uberon (`UBERON:`) ontology is used. |
| `obs["disease"]` | `string` | (*Optional*) Information on any disease or pathological condition associated with the cell or donor. |
| `obs["disease_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the disease, enabling standardized disease classification and referencing. Must be a term from the Mondo Disease Ontology (`MONDO:`) ontology term, or `PATO:0000461` from the Phenotype And Trait Ontology (`PATO:`). |
| `obs["donor_id"]` | `string` | (*Optional*) Identifier for the donor from whom the cell sample is obtained. |
| `obs["is_primary_data"]` | `boolean` | (*Optional*) Indicates whether the data is primary (directly obtained from experiments) or has been computationally derived from other primary data. |
| `obs["organism"]` | `string` | (*Optional*) Organism from which the cell sample is obtained. |
| `obs["organism_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the organism, providing a standardized reference for the organism. Must be a term from the NCBI Taxonomy Ontology (`NCBITaxon:`) which is a child of `NCBITaxon:33208`. |
| `obs["self_reported_ethnicity"]` | `string` | (*Optional*) Ethnicity of the donor as self-reported, relevant for studies considering genetic diversity and population-specific traits. |
| `obs["self_reported_ethnicity_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the self-reported ethnicity, providing a standardized reference for ethnic classifications. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Ancestry Ontology (`HANCESTRO:`) is used. |
| `obs["sex"]` | `string` | (*Optional*) Biological sex of the donor or source organism, crucial for studies involving sex-specific traits or conditions. |
| `obs["sex_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the biological sex, ensuring standardized classification of sex. Only `PATO:0000383`, `PATO:0000384` and `PATO:0001340` are allowed. |
| `obs["suspension_type"]` | `string` | (*Optional*) Type of suspension or medium in which the cells were stored or processed, important for understanding cell handling and conditions. |
| `obs["tissue"]` | `string` | (*Optional*) Specific tissue from which the cells were derived, key for context and specificity in cell studies. |
| `obs["tissue_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the tissue, providing a standardized reference for the tissue type. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`. |
| `obs["tissue_general"]` | `string` | (*Optional*) General category or classification of the tissue, useful for broader grouping and comparison of cell data. |
| `obs["tissue_general_ontology_term_id"]` | `string` | (*Optional*) Ontology term identifier for the general tissue category, aiding in standardizing and grouping tissue types. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`. |
| `obs["batch"]` | `string` | (*Optional*) A batch identifier. This label is very context-dependent and may be a combination of the tissue, assay, donor, etc. |
| `obs["soma_joinid"]` | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the cell. |
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature, usually a ENSEMBL gene id. |
| `var["feature_name"]` | `string` | A human-readable name for the feature, usually a gene symbol. |
| `var["soma_joinid"]` | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the feature. |
| `var["hvg"]` | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `var["hvg_score"]` | `double` | A score for the feature indicating how highly variable it is. |
| `obsm["X_pca"]` | `double` | The resulting PCA embedding. |
| `obsp["knn_distances"]` | `double` | K nearest neighbors distance matrix. |
| `obsp["knn_connectivities"]` | `double` | K nearest neighbors connectivities matrix. |
| `varm["pca_loadings"]` | `double` | The PCA loadings matrix. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["normalized"]` | `integer` | Normalized expression values. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. This is different from the `obs.dataset_id` field, which is the identifier for the dataset from which the cell data is derived. |
| `uns["dataset_name"]` | `string` | A human-readable name for the dataset. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

