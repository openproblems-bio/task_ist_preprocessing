# Load required packages
library(data.table)
library(yaml)
library(anndata)
library(FRmatch)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(jsonlite)
library(here)

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
    filter_size = NULL,
    filter_fscore = NULL,
    filter_nomarker = NULL,
    add_pseudo_marker = NULL,
    pseudo_expr = NULL,
    subsamp_size = NULL,
    subsamp_iter = NULL,
    subsamp_seed = NULL,
    numCores = NULL,
    per_gene_correction = NULL,
    per_gene_layer = NULL,
    cell_type_column = NULL,
    leiden_res = NULL,
    n_pcs = NULL,
    cell_id = NULL
  )
  for (i in seq_along(args)) {
    if (args[i] %in% c("-s", "--spatial")) {
      arg_dict$spatial <- args[i + 1]
    } else if (args[i] %in% c("-d", "--dissociated")) {
      arg_dict$dissociated <- args[i + 1]
    } else if (args[i] %in% c("-o", "--output")) {
      arg_dict$output <- args[i + 1]
    }  else if (args[i] %in% c("-p", "--hyperparameters")) {
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
    
    
    if ("filter_size" %in% names(params)) {
      arg_dict$filter_size <- as.integer(params$filter_size)
    }
    if ("filter_fscore" %in% names(params)) {
      arg_dict$filter_fscore <- as.numeric(params$filter_fscore)
    }
    if ("filter_nomarker" %in% names(params)) {
      arg_dict$filter_nomarker <- as.logical(params$filter_nomarker)
    }
    if ("add_pseudo_marker" %in% names(params)) {
      arg_dict$add_pseudo_marker <- as.logical(params$add_pseudo_marker)
    }
    if ("pseudo_expr" %in% names(params)) {
      arg_dict$pseudo_expr <- as.integer(params$pseudo_expr)
    }
    if ("subsamp_size" %in% names(params)) {
      arg_dict$subsamp_size <- as.integer(params$subsamp_size)
    }
    if ("subsamp_iter" %in% names(params)) {
      arg_dict$subsamp_iter <- as.integer(params$subsamp_iter)
    }
    if ("subsamp_seed" %in% names(params)) {
      arg_dict$subsamp_seed <- as.integer(params$subsamp_seed)
    }
    if ("numCores" %in% names(params)) {
      arg_dict$numCores <- as.integer(params$numCores)
    }
    if ("cell_type" %in% names(params)) {
      arg_dict$cell_type_column <- params$cell_type
    }
    if ("leiden_res" %in% names(params)) {
      arg_dict$leiden_res <- as.numeric(params$leiden_res)
    }
    if ("n_pcs" %in% names(params)) {
      arg_dict$n_pcs <- as.integer(params$n_pcs)
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

  hparams_defaults <- defaults$frmatch 
  gparams_defaults <- defaults$annotation_params
  
  # Parse hyperparameters
  hyperparams <- list(
  feature_selection = "query.genes", 
  filter_size = if (!is.null(args$filter_size)) as.integer(args$filter_size) else hparams_defaults$filter_size,
  filter_fscore = if (!is.null(args$filter_fscore)) as.numeric(args$filter_fscore) else hparams_defaults$filter_fscore,
  filter_nomarker = if (!is.null(args$filter_nomarker)) as.logical(args$filter_nomarker) else hparams_defaults$filter_nomarker,
  add_pseudo_marker = if (!is.null(args$add_pseudo_marker)) as.logical(args$add_pseudo_marker) else hparams_defaults$add_pseudo_marker,
  pseudo_expr = if (!is.null(args$pseudo_expr)) as.integer(args$pseudo_expr) else hparams_defaults$pseudo_expr,
  subsamp_size = if (!is.null(args$subsamp_size)) as.integer(args$subsamp_size) else hparams_defaults$subsamp_size,
  subsamp_iter = if (!is.null(args$subsamp_iter)) as.integer(args$subsamp_iter) else hparams_defaults$subsamp_iter,
  subsamp_seed = if (!is.null(args$subsamp_seed)) as.integer(args$subsamp_seed) else hparams_defaults$subsamp_seed,
  numCores = if (!is.null(args$numCores)) as.integer(args$numCores) else hparams_defaults$numCores,
  prefix_q = "",
  prefix_ref = "",
  cell_type_column = if (!is.null(args$cell_type_column)) args$cell_type_column else hparams_defaults$cell_type_column,
  leiden_res = if (!is.null(args$leiden_res)) as.numeric(args$leiden_res) else hparams_defaults$leiden_res,
  n_pcs = if (!is.null(args$n_pcs)) as.integer(args$n_pcs) else hparams_defaults$n_pcs,
  cell_id = if (!is.null(args$cell_id)) args$cell_id else hparams_defaults$cell_id
  )


  groupparams <- list(
    per_gene_correction = if (!is.null(args$per_gene_correction)) args$per_gene_correction else gparams_defaults$per_gene_correction,
    per_gene_layer = if (!is.null(args$per_gene_layer)) args$per_gene_layer else gparams_defaults$per_gene_layer
    
  )
  # Print parsed parameters
  cat("Hyperparameters:", paste(hyperparams, collapse=", "), "\n")
  cat("Group parameters:", paste(groupparams, collapse=", "), "\n")
  
  # Read data
  adata_sc <- read_h5ad(args$dissociated)
  adata <- read_h5ad(args$spatial)
  
  # make data objects
  

  # single cell data object
  sc_df = as.data.frame(as.matrix(adata_sc$X))
  sc_df <- cbind(rownames(adata_sc$obs), sc_df)
  colnames(sc_df)[1] <- "Sample"


  cluster_labels <- as.character(adata_sc$obs[[hyperparams$cell_type_column]])

  # Create the data frame
  cell_cluster_labels_sc <- data.frame(
    Sample = rownames(adata_sc$obs),
    Cluster = cluster_labels
  )

  object_sc <- make_data_object(dat = sc_df, 
                            tab = cell_cluster_labels_sc,
                            markers = colnames(adata$X))  # markers are genes from spatial data
          
  #spatial data object    
  # run clustering first
  print('Initial Clustering')
  adata_Seurat <- CreateSeuratObject(counts = tcrossprod(adata$X), project = "spatial", min.cells = 3, min.features = 200)
  
  

  all.genes <- rownames(adata_Seurat)
  adata_Seurat <- ScaleData(adata_Seurat,features = all.genes)

  adata_Seurat <- RunPCA(adata_Seurat, features = all.genes)
  adata_Seurat <- FindNeighbors(adata_Seurat, dims = 1:hyperparams$n_pcs)
  adata_Seurat <- FindClusters(adata_Seurat, resolution = hyperparams$leiden_res, algorithm = 4)
  
  
  cell_cluster_labels_spatial <- data.frame(Sample = adata$obs[[hyperparams$cell_id]], Cluster = Idents(adata_Seurat))


  #create object
  sp_df = as.data.frame(as.matrix(adata$X))
  sp_df <- cbind(adata$obs[[hyperparams$cell_id]], sp_df)
  colnames(sp_df)[1] <- "Sample"
  object_spatial <- make_data_object(dat = sp_df,
                            tab = cell_cluster_labels_spatial,
                            markers = colnames(adata$X)) 

  # Run cell type annotation
  print('Mapping cell types')
  result <- FRmatch_cell2cluster(object_spatial, object_sc,
                      feature.selection=hyperparams$feature_selection, 
                      filter.size=hyperparams$filter_size, filter.fscore=hyperparams$filter_fscore, filter.nomarker=hyperparams$filter_nomarker, #filtering clusters
                      add.pseudo.marker=hyperparams$add_pseudo_marker, pseudo.expr=hyperparams$pseudo_expr, #adding pseudo marker
                      subsamp.size=hyperparams$subsamp_size, subsamp.iter=hyperparams$subsamp_iter, subsamp.seed=hyperparams$subsamp_seed, #subsampling
                      numCores=hyperparams$numCores, prefix=c(hyperparams$prefix_q, hyperparams$prefix_ref))
  
  # check result type
  annotation_df <- result$cell2cluster
  #rename columns
  colnames(annotation_df)[colnames(annotation_df) == "match"] <- "celltype"
  colnames(annotation_df)[colnames(annotation_df) == "query.cell"] <- "cell_id"
 
  # Keep only 'cell_id', 'celltype', and 'score' columns
  annotation_df <- annotation_df[, c('cell_id', 'celltype', 'score')]
  # Save annotation
  write.csv(annotation_df, file=args$output, row.names = FALSE)
}




# Entry point of the script
main <- function() {
  
  print("Starting script execution...")
  args <- parse_args()
  print("Arguments parsed successfully.")
  annotate_cells(args)
  print("Script execution completed.")
}

# Call main function
main()
