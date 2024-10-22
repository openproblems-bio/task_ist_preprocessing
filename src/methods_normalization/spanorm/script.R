#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("remotes")
#BiocManager::install("bhuvad/SpaNorm")
#BiocManager::install("zellkonverter")
library(SpaNorm)
library(SingleCellExperiment)
library(SpatialExperiment)
library(zellkonverter)

par <- list(
  "input" = 'spatial_aggregated_counts.h5ad',
  "output" = 'spatial_spanormed_counts.h5ad'
)

sce <- readH5AD(par$input)
sce <- as(sce, "SpatialExperiment")

centroid_x <- colData(sce)$centroid_x
centroid_y <- colData(sce)$centroid_y

spatial_coords <- matrix(c(centroid_x, centroid_y), ncol = 2)

spatialCoords(sce) <- spatial_coords

keep = filterGenes(sce, 0.05)
sce <- sce[keep, ]
result <- SpaNorm(sce)

main_assay <- assay(result, "X")
row_data <- rowData(result)
col_data <- colData(result)

final_sce <- SingleCellExperiment(
  assays = list(X = main_assay),
  rowData = row_data,
  colData = col_data
)

zellkonverter::writeH5AD(final_sce, par$output)

