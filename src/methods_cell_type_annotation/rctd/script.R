library(spacexr)
library(Matrix)
library(SingleCellExperiment)
# library(SpatialExperiment)
library(anndataR)

## VIASH START
par <- list(
  "input_spatial_normalized_counts" = "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/spatial_aggregated_counts.h5ad",
  "input_scrnaseq_reference"= "task_ist_preprocessing/resources_test/task_ist_preprocessing/mouse_brain_combined/scrnaseq_reference.h5ad",
  "output" = "task_ist_preprocessing/tmp/spatial_types.h5ad",
  "gene_cutoff" = 0.0,
  "fc_cutoff" = 0.1,
  "gene_cutoff_reg" = 0.0,
  "fc_cutoff_reg" = 0.1,
  "umi_min" = 20,
  "umi_min_sigma" = 20
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

# --- Diagnostics: confirm what RCTD actually reads. The reference itself has
#     ~475/477 panel genes with clear cross-cell-type fold-change, so a "fewer
#     than 10 DE genes" failure means the data reaching RCTD is degenerate:
#     cell types collapsed to ~1, non-integer/normalized counts, or few genes
#     shared with the spatial puck. This block surfaces which. ---
cat("=== RCTD input diagnostics ===\n")
cat(sprintf("spatial puck: %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("reference (>=25 cells/type): %d genes x %d cells, %d cell types\n",
            nrow(ref_counts), ncol(ref_counts), length(unique(as.character(cell_types)))))
print(utils::head(sort(table(as.character(cell_types)), decreasing = TRUE), 8))
cat(sprintf("genes shared reference<->spatial: %d\n",
            length(intersect(rownames(ref_counts), rownames(counts)))))
.rcv <- if (methods::is(ref_counts, "sparseMatrix")) ref_counts@x else as.numeric(ref_counts)
if (length(.rcv)) cat(sprintf("reference counts: min=%.4g max=%.4g mean=%.4g all-integer=%s\n",
            min(.rcv), max(.rcv), mean(.rcv), isTRUE(all.equal(.rcv, round(.rcv)))))
cat("==============================\n")

# spacexr::check_cell_types() rejects cell-type names containing '/'
# (e.g. "Ciliated/secretory cells", "T/NK lineage"). Sanitize the factor levels
# before building the Reference, and keep a map (safe name -> original name) so
# the predicted labels can be restored below to match the scRNAseq reference.
cell_type_name_map <- levels(cell_types)
names(cell_type_name_map) <- gsub("/", "_", levels(cell_types))
if (any(duplicated(names(cell_type_name_map)))) {
  names(cell_type_name_map) <- make.unique(names(cell_type_name_map))
}
levels(cell_types) <- names(cell_type_name_map)

reference <- Reference(ref_counts, cell_types, min_UMI = 0)

# check cores
cores <- 1
if ("cpus" %in% names(meta) && !is.null(meta$cpus)) cores <- meta$cpus
cat(sprintf("Number of cores: %s\n", cores))

# Run the algorithm
# NOTE: RCTD's default DE-gene / UMI thresholds are tuned for whole-transcriptome
# references and produce "fewer than 10 regression differentially expressed genes"
# on small iST panels. The relaxed thresholds below are passed from the config.
myRCTD <- create.RCTD(
  puck, reference, max_cores = cores,
  gene_cutoff = par$gene_cutoff, fc_cutoff = par$fc_cutoff,
  gene_cutoff_reg = par$gene_cutoff_reg, fc_cutoff_reg = par$fc_cutoff_reg,
  UMI_min = par$umi_min, UMI_min_sigma = par$umi_min_sigma
)
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
# Restore the original cell-type names (undo the '/' sanitization). "None_sp"
# and anything not in the map are left unchanged.
predicted <- as.character(spatial_cell_types)
mapped <- unname(cell_type_name_map[predicted])
predicted[!is.na(mapped)] <- mapped[!is.na(mapped)]
colData(sce)[names(spatial_cell_types),"cell_type"] <- predicted

# Write the final object to h5ad format
# set to 'w', is this ok?
dir.create(dirname(par$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(sce, par$output, mode = "w")
