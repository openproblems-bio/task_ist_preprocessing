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
#I put random coordinates 
num_cells <- ncol(sce)
set.seed(42)
random_x <- runif(num_cells, min = 0, max = 100) 
random_y <- runif(num_cells, min = 0, max = 100)  
random_coords <- matrix(c(random_x, random_y), ncol = 2)
spatialCoords(sce) <- random_coords

keep = filterGenes(sce, 0.05)
sce <- sce[keep, ]
result <- SpaNorm(sce)
zellkonverter::writeH5AD(as(result, 'SingleCellExperiment'), par$output)

