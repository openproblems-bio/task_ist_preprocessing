
# Load required packages
library(data.table)
library(yaml)
library(anndata)
library(jsonlite)
library(here)
suppressPackageStartupMessages({
  library(scrattch.mapping)
  library(scrattch.taxonomy)   
})

print("Libraries loaded.")


# Define argument parser function
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6) {
    stop("Error: Missing required arguments. Usage: Rscript your_script.R -s spatial_file -d dissociated_file -o output_file [-p hyperparameters -g group_parameters]")
  }
  p=NULL
  g= NULL
  arg_dict <- list(
    spatial = NULL,
    dissociated = NULL,
    output = NULL,
    temp_dir=NULL,
    method = NULL, 
    cell_type_col = NULL,
    cell_id= NULL,
    per_gene_correction = NULL,
    per_gene_layer = NULL
    
   
   
  )
  for (i in seq_along(args)) {
    if (args[i] %in% c("-s", "--spatial")) {
      arg_dict$spatial <- args[i + 1]
    } else if (args[i] %in% c("-d", "--dissociated")) {
      arg_dict$dissociated <- args[i + 1]
    }  else if (args[i] %in% c("-o", "--output")) {
      arg_dict$output <- args[i + 1]
    }  else if (args[i] %in% c("-t", "-- temporary output dir")) {
      arg_dict$temp_dir <- args[i + 1]
    } else if (args[i] %in% c("-p", "--hyperparamters")) {
      p <- args[i + 1]
    }  else if (args[i] %in% c("-g", "--group-paramters")) {
      g <- args[i + 1]
    }
  }

  # Check if -p parameter exists
  if (!is.null(p)) {
    # Parse JSON string
    json_string <- gsub("\\'", '"', p )  # Replace single quotes with double quotes
    params <- fromJSON(json_string)
    
    if ("cell_type_col" %in% names(params)) {
      arg_dict$cell_type_col <- params$cell_type_col
    }
    if ("method" %in% names(params)) {
      arg_dict$method <- params$method
    }
    if ("cell_id" %in% names(params)) {
      arg_dict$cell_id <- params$cell_id
    }
    
  }
  if (!is.null(g)){
    # Parse JSON string
    json_string <- gsub("\\'", '"', g )  # Replace single quotes with double quotes
    params <- fromJSON(json_string)
    if ("per_gene_correction" %in% names(params)) {
      arg_dict$per_gene_correction <- as.logical(params$per_gene_correction)
    }
    if ("per_gene_layer" %in% names(params)) {
      arg_dict$per_gene_layer <- params$per_gene_layer
    }

  }

  return(arg_dict)
}


# Function to create output directory if it doesn't exist
create_output_dir <- function(output_path) {
  dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)
}


# Main function to run annotation
annotate_cells <- function(args) {
  create_output_dir(dirname(args$output))
  
  # Load default hyperparameters from YAML config
  #defaults <- read_yaml("configs/defaults.yaml")
  defaults <- read_yaml(file.path(here("configs"), "defaults.yaml"))
  hparams_defaults <- defaults$scrattchmapping
  gparams_defaults <- defaults$annotation_params
  
  # Parse hyperparameters
  hyperparams <- list(
    cell_type_col = if (!is.null(args$cell_type_col)) args$cell_type_col else hparams_defaults$cell_type_col,
    method = if (!is.null(args$method)) args$method else hparams_defaults$method,
    cell_id = if (!is.null(args$cell_id)) args$cell_id else hparams_defaults$cell_id
    
  )
  groupparams <- list(
    per_gene_correction = if (!is.null(args$per_gene_correction)) args$per_gene_correction else gparams_defaults$per_gene_correction,
    per_gene_layer = if (!is.null(args$per_gene_layer)) args$per_gene_layer else gparams_defaults$per_gene_layer
    
  )
  
  # Print parsed parameters
  cat("Hyperparameters:", paste(hyperparams, collapse=", "), "\n")
  cat("Group parameters:", paste(groupparams, collapse=", "), "\n")

  if (hyperparams$method == "correlation") {
      corr.map <- TRUE
      tree.map <- TRUE
      seurat.map <- FALSE
  } else if (hyperparams$method == "tree") {
      corr.map <- FALSE
      tree.map <- TRUE
      seurat.map <- FALSE
  } else if (hyperparams$method == "seurat") {
      corr.map <- FALSE
      tree.map <- TRUE
      seurat.map <- TRUE
  } else {
      # Handle the case if hyperparams$method is none of the specified values
      stop("Invalid method specified.")
  }

  


  # Read data
  adata_sc <- read_h5ad(args$dissociated)
  adata <- read_h5ad(args$spatial)
  spatial_data = t(as.matrix(adata$X)) # transpose matrix : genes as rows and cells as columns

  print('Deleting')
  # delete sc celltype with less than 2 cells
  category_counts <- table(adata_sc$obs[[hyperparams$cell_type_col]])
  categories_to_keep <- names(category_counts[category_counts >= 2])
  rows_to_delete <- which(!adata_sc$obs[[hyperparams$cell_type_col]] %in% categories_to_keep)
  taxonomy.anno <- adata_sc$obs[-rows_to_delete, ]
  taxonomy.counts <- t(as.matrix(adata_sc$X[-rows_to_delete, ])) # Transpose matrix: genes as rows and cells as columns
  print('done')



  ## Ensure 'cluster' field exists, as required by scrattch.taxonomy.
  taxonomy.anno$cluster = taxonomy.anno[[hyperparams$cell_type_col]]

  ## Compute top 1000 binary marker genes for clusters (or use a pre-existing vector)
  #binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster, 1000)
  binary.genes <- intersect(rownames(spatial_data),rownames(taxonomy.counts))
  ## Compute UMAP coordinates (or use precomputed coordinates)
  # print('Computing PCs')
  # pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
  # print('Computing umap')
  # umap.coords = umap(pcs[,1:30])$layout

  ## Set rownames to your annotation and UMAP data.frames with sample identifiers (Required!)
  ##rownames(taxonomy.anno) = taxonomy.anno$sample_name
  umap.coords <- data.frame(
  x = taxonomy.anno$x,  
  y = taxonomy.anno$y   
  )
  rownames(umap.coords) = colnames(taxonomy.counts)
  taxonomyDir = args$temp_dir

  ## Build Shiny taxonomy 
  print('Creating object')
  AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
                  meta.data = taxonomy.anno,
                  feature.set = binary.genes,
                  umap.coords = umap.coords,
                  taxonomyName = "sc_data", ## NEW!
                  taxonomyDir = taxonomyDir,
                  subsample=2000)

  ## Add markers to dendrogram
  write.csv(AIT.anndata$obs,file = 'temp.csv')
  print('Adding dendogram markers')
  AIT.anndata = addDendrogramMarkers(AIT.anndata = AIT.anndata)


  print("Mapping cell types...")
  mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata,
                                query.data=spatial_data,
                                label.cols="cluster_label", ## Which obs in AIT.anndata contain annotations to map. E.g. "class", "subclass", etc.
                                corr.map=corr.map,
                                tree.map=tree.map,
                                seurat.map=seurat.map)

  ## Extract mapping results from S4 mappingClass
  annotation_df = getMappingResults(mapping.anno)

  ## Extract tree mapping bootstraping table (We will improve this in the near future.)
  #tree.bootstraps = mapping.anno@detailed_results[["tree"]]


  annotation_df$cell_id <- adata$obs[[hyperparams$cell_id]]
  # Reorder columns in annotation_df with cell_id as the first column
  annotation_df <- annotation_df[, c("cell_id", setdiff(names(annotation_df), "cell_id"))]
  #rename columns
  
  
  if (hyperparams$method == "correlation") {
    colnames(annotation_df)[colnames(annotation_df) == "cluster_Corr"] <- "celltype"
    colnames(annotation_df)[colnames(annotation_df) == "score.Corr"] <- "score"
  } else if (hyperparams$method == "tree") {
    colnames(annotation_df)[colnames(annotation_df) == "cluster_Tree"] <- "celltype"
    colnames(annotation_df)[colnames(annotation_df) == "score.Tree"] <- "score"
  } else if (hyperparams$method == "seurat") {
    colnames(annotation_df)[colnames(annotation_df) == "cluster_Seurat"] <- "celltype"
    colnames(annotation_df)[colnames(annotation_df) == "score.Seurat"] <- "score"
  }
  

  # Keep only 'cell_id', 'celltype', and 'score' columns
  annotation_df <- annotation_df[, c('cell_id', 'celltype', 'score')]

  # Save annotation
  write.csv(annotation_df, file=args$output, row.names = FALSE)


}




# Entry point of the script
main <- function() {
  args <- parse_args()
  annotate_cells(args)
}

# Call main function
main()