library(spacexr)
library(Matrix)
library(SingleCellExperiment)
library(anndataR)

## VIASH START
par <- list(
  "input_spatial_normalized_counts" = "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad",
  "input_scrnaseq_reference"= "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad",
  "output" = "task_ist_preprocessing/tmp/spatial_types.h5ad"
)

meta <- list(
  'cpus': 4,
)

## VIASH END

# Read the input h5ad file and convert to SingleCellExperiment
sce <- read_h5ad(par$input_spatial_normalized_counts, as = "SingleCellExperiment")

# Extract spatial coordinates and counts matrix
centroid_x <- colData(sce)$centroid_x
centroid_y <- colData(sce)$centroid_y
coords <- data.frame(centroid_x, centroid_y)
counts <- assay(sce,"counts")
rownames(coords) <- colData(sce)$cell_id
puck <- SpatialRNA(coords, counts)

# Read reference scrnaseq
ref <- read_h5ad(par$input_scrnaseq_reference, as = "SingleCellExperiment")

#filter reference cell types to those with >25 cells
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
myRCTD <- create.RCTD(puck, reference, max_cores = cores)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# Extract results
results <- myRCTD@results
spatial_cell_types <- results$results_df$first_type 
# Include None Spatial cell type for the "reject" cells
levels(spatial_cell_types) <- c(levels(spatial_cell_types), "None_sp")
spatial_cell_types[results$results_df$spot_class == "reject"] <- "None_sp"
names(spatial_cell_types) <- rownames(results$results_df)

#
colData(sce)$cell_type <- "None_sp"
colData(sce)[names(spatial_cell_types),"cell_type"] <- as.character(spatial_cell_types)

# Write the final object to h5ad format
dir.create(dirname(par$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(sce, par$output)
