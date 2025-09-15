library(SpaNorm)
library(SingleCellExperiment)
library(SpatialExperiment)
library(zellkonverter)

## VIASH START
par <- list(
  "input_spatial_aggregated_counts" = 'resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad',
  "output" = 'tmp/spatial_spanormed_counts.h5ad'
)
## VIASH END

# Read the input h5ad file and convert to SingleCellExperiment
sce <- readH5AD(par$input_spatial_aggregated_counts)
# Convert to SpatialExperiment for SpaNorm
sce <- as(sce, "SpatialExperiment")

# Extract spatial coordinates
centroid_x <- colData(sce)$centroid_x
centroid_y <- colData(sce)$centroid_y

# Create spatial coordinates matrix for SpaNorm
spatial_coords <- matrix(c(centroid_x, centroid_y), ncol = 2)

# Set spatial coordinates in the SpatialExperiment object
spatialCoords(sce) <- spatial_coords

# Apply SpaNorm normalization to the spatial data
set.seed(777)

# Calculate sample.p based on number of cells (following recommendations of https://github.com/bhuvad/SpaNorm/issues/10#issuecomment-2918248535)
n_cells <- ncol(sce)
if (n_cells > 100000) {
  sample_p <- 25000 / n_cells
} else {
  sample_p <- 0.25
}

result <- SpaNorm(sce, sample.p = sample_p)

# Get the normalized matrix from SpaNorm result (log-transformed normalized counts)
normalized_matrix <- assay(result, "logcounts")

# Create final SCE with all original layers preserved
final_sce <- SingleCellExperiment(
  assays = assays(sce),  # Preserve all original assays
  rowData = rowData(sce),
  colData = colData(sce)
)

# Add the normalized matrix as a new layer called 'normalized'
assay(final_sce, "normalized") <- normalized_matrix

# Write the final object to h5ad format
zellkonverter::writeH5AD(final_sce, par$output)
