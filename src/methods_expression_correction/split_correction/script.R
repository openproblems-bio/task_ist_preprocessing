library(spacexr)
library(Matrix)
library(SingleCellExperiment)
library(anndataR)
library(SPLIT)
library(Seurat)
library(scuttle)

## VIASH START
par <- list(
  "input_spatial_with_cell_types" = "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad",
  "input_scrnaseq_reference"= "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad",
  "output" = "task_ist_preprocessing/tmp/split_corrected.h5ad",
  "keep_all_cells" = FALSE,
)

meta <- list(
  'cpus': 4,
)

## VIASH END

# Read the input h5ad file and convert to SingleCellExperiment and Seurat
sce <- read_h5ad(par$input_spatial_with_cell_types, as = "SingleCellExperiment")
xe <- read_h5ad(par$input_spatial_with_cell_types, as = "Seurat")

# filter out 0 cells
if (!par$keep_all_cells) {
  cat("Filtering cells with 0 counts\n")
  sce <- sce[, colSums(counts(sce)) > 0]
  xe <- subset(xe, subset = nCount_RNA > 0)
}

# Extract spatial coordinates and counts matrix
centroid_x <- colData(sce)$centroid_x
centroid_y <- colData(sce)$centroid_y
coords <- data.frame(centroid_x, centroid_y)
counts <- assay(sce, "counts")
rownames(coords) <- colData(sce)$cell_id
puck <- SpatialRNA(coords, counts)

# Read reference scrnaseq
ref <- read_h5ad(par$input_scrnaseq_reference, as = "SingleCellExperiment")

#filter reference cell types to those with >25 cells (minimum for RCTD)
valid_celltypes <- names(table(colData(ref)$cell_type))[table(colData(ref)$cell_type) >= 25] 
filtered_ref <- ref[,colData(ref)$cell_type %in% valid_celltypes]

ref_counts <- assay(filtered_ref, "counts")
# factor to drop filtered cell types
colData(filtered_ref)$cell_type <- factor(colData(filtered_ref)$cell_type)
cell_types <- colData(filtered_ref)$cell_type
names(cell_types) <- colnames(ref_counts)
reference <- Reference(ref_counts, cell_types, min_UMI = 0)

# check cores
cores <- 1
if ("cpus" %in% names(meta) && !is.null(meta$cpus)) cores <- meta$cpus
cat(sprintf("Number of cores: %s\n", cores))

# Run the algorithm
cat("Running RCTD\n")
myRCTD <- create.RCTD(puck, reference, max_cores = cores)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# Get the "spot_class" annotation from RCTD
# cat("Saving RCTD spot_class\n")
# results <- myRCTD@results
# rctd_spot_class <- results$results_df$spot_class
# names(rctd_spot_class) <- rownames(results$results_df)
# colData(sce)$RCTD_class <- "not_included"
# colData(sce)[names(rctd_spot_class),"RCTD_class"] <- as.character(rctd_spot_class)

# Post-process RCTD output
RCTD <- SPLIT::run_post_process_RCTD(myRCTD)

# Run SPLIT purification
cat("Running SPLIT\n")
res_split <- SPLIT::purify(
  counts = GetAssayData(xe, assay = 'RNA', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = TRUE # optional
)


# create corrected counts layer in original SingleCell object
cat("Normalizing counts\n")

# First copy in counts
assay(sce, "corrected_counts") <- assay(sce, "counts")

# Then, replace only the updated cells
assay(sce, "corrected_counts")[rownames(res_split$purified_counts), colnames(res_split$purified_counts)] <- res_split$purified_counts

# Library size normalization - see note in resolVI
size_factors <- librarySizeFactors(assay(sce, "corrected_counts"))
assay(sce, "normalized") <- assay(logNormCounts(sce, size_factors=size_factors, assay.type = "corrected_counts"),"logcounts")

# Write the final object to h5ad format
cat("Writing to h5ad\n")
dir.create(dirname(par$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(sce, par$output, mode = "w")
