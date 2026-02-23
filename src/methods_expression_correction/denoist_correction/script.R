library(Matrix)
library(DenoIST)
library(SpatialExperiment)
library(SingleCellExperiment)
library(anndataR)
library(scuttle)

## VIASH START
par <- list(
  "input_spatial_with_cell_types" = "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad",
  "input_tx" = "mouse_combined_transcripts.csv",
  "output" = "task_ist_preprocessing/tmp/split_corrected.h5ad",
#   "keep_all_cells" = FALSE,
  "distance" = 50,
  "nbins" = 200,
)

meta <- list(
  'cpus': 4,
)

## VIASH END

# Read the input h5ad file and convert to SingleCellExperiment and Seurat
sce <- read_h5ad(par$input_spatial_with_cell_types, as = "SingleCellExperiment")
spe <- SpatialExperiment(
    assay = counts(sce),
    colData = colData(sce),
    spatialCoordsNames = c("centroid_x", "centroid_y"))

tx <- read.csv(par$input_tx)

# filter out 0 cells
# if (!par$keep_all_cells) {
#   cat("Filtering cells with 0 counts\n")
#   sce <- sce[, colSums(counts(sce)) > 0]
#   xe <- subset(xe, subset = nCount_RNA > 0)
# }


# check cores
cores <- 1
if ("cpus" %in% names(meta) && !is.null(meta$cpus)) cores <- meta$cpus
cat(sprintf("Number of cores: %s\n", cores))

# Run the algorithm

res <- denoist(mat = spe,
              tx = tx,
              feature_label = "feature_name",
              coords = NULL,
              distance = par$distance, nbins = par$nbins, cl = cores) #TODO add in params

# format name
corrected_counts <- res$adjusted_counts

# create corrected counts layer in original SingleCell object
cat("Normalizing counts\n")

# First copy in counts
assay(sce, "corrected_counts") <- assay(sce, "counts")

# Then, replace only the updated cells
assay(sce, "corrected_counts")[rownames(corrected_counts), colnames(corrected_counts)] <- corrected_counts

# Library size normalization - see note in resolVI
size_factors <- librarySizeFactors(assay(sce, "corrected_counts"))
assay(sce, "normalized") <- assay(logNormCounts(sce, size_factors=size_factors, assay.type = "corrected_counts"),"logcounts")

# Write the final object to h5ad format
cat("Writing to h5ad\n")
dir.create(dirname(par$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(sce, par$output, mode = "w")