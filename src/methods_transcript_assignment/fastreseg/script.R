devtools::install_github("Nanostring-Biostats/FastReseg", build_vignettes = FALSE, ref = "main")

library(FastReseg)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  stop("No arguments provided")
}

# Access individual arguments
path_to_counts <- args[1]
path_to_transcripts <- args[2]
path_to_cell_annot <- args[3]

path_to_cell_ids_out  <- args[4]
path_to_gene_names_out <- args[5] 
path_to_transcripts_out <- args[6]


### reading in data
count_df <- as.matrix(read.csv(path_to_counts, row.names = 1))
head(count_df)

cell_df <- read.csv(path_to_cell_annot, row.names = 1)
cell_annot <- cell_df[row.names(count_df),]
cell_annot <- setNames(cell_annot, row.names(count_df))

transcriptDF <- read.csv(path_to_transcripts)
head(transcriptDF)


##run pipeline
refineAll_res_one_FOC <- fastReseg_full_pipeline(
  counts = count_df,
  clust = cell_annot,
  
  # Similar to `runPreprocess()`, one can use `clust = NULL` if providing `refProfiles`
  
  transcript_df = transcriptDF, #
  transDF_fileInfo = NULL,
  pixel_size = 0.18,
  zstep_size = 0.8,
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = "UMI_cellID",
  spatLocs_colns = c("x","y","z"),
  extracellular_cellID = c(-1),
  
  # Similar to `runPreprocess()`, one can set various cutoffs to NULL for automatic calculation from input data
  
  # distance cutoff for neighborhood searching at molecular and cellular levels, respectively
  molecular_distance_cutoff = 2.7, 
  cellular_distance_cutoff = NULL, 
  
  # cutoffs for transcript scores and number for cells under each cell type
  score_baseline = NULL,
  lowerCutoff_transNum = NULL,
  higherCutoff_transNum= NULL,
  imputeFlag_missingCTs = TRUE,
  
  # Settings for error detection and correction, refer to `runSegRefinement()` for more details
  flagCell_lrtest_cutoff = 5, # cutoff to flag for cells with strong spatial dependcy in transcript score profiles
  svmClass_score_cutoff = -2,   # cutoff of transcript score to separate between high and low score classes
  groupTranscripts_method = "dbscan",
  spatialMergeCheck_method = "leidenCut", 
  cutoff_spatialMerge = 0.5, # spatial constraint cutoff for a valid merge event
  
  path_to_output = "res2_multiFiles",
  save_intermediates = TRUE, # flag to return and write intermediate results to disk
  return_perCellData = TRUE, # flag to return per cell level outputs from updated segmentation 
  combine_extra = FALSE # flag to include trimmed and extracellular transcripts in the exported `updated_transDF.csv` files 
)

if(is.null(refineAll_res_one_FOC$updated_transDF_list[[1]])){
  print("no transcripts assigned after refinement")
  cell_df$UMI_cellID <- as.integer(row.names(cell_df))
  transcriptDF <- left_join(cell_df, transcriptDF)
  names(transcriptDF)[names(transcriptDF) == "UMI_cellID"] <- "updated_cellID"
  names(transcriptDF)[names(transcriptDF) == "ct_ssam"] <- "updated_celltype"
  write.csv(transcriptDF, file = path_to_transcripts_out)
} else {
  write.csv(refineAll_res_one_FOC$updated_transDF_list[[1]], file = path_to_transcripts_out)
}
### export outputs
write.csv(row.names(refineAll_res_one_FOC$updated_perCellExprs), file = path_to_gene_names_out)
write.csv(refineAll_res_one_FOC$updated_perCellDT, file = path_to_cell_ids_out)
