# Load required packages
library(data.table)
library(yaml)
library(anndata)
library(jsonlite)
library(here)
suppressPackageStartupMessages({
  library(mfishtools)    
  library(matrixStats)   
})
options(stringsAsFactors = FALSE)  
print("Libraries loaded.")

cellToClusterMapping_byCor_ <- function(medianDat,
                                       mapDat,
                                       refDat = NA,
                                       clusters = NA,
                                       genesToMap = rownames(mapDat),
                                       use = "p",
                                       method = "p",
				                               returnCor=TRUE,
                                       ...) {
  corVar <- corTreeMapping(
    medianDat = medianDat,
    mapDat = mapDat, refDat = refDat, clusters = clusters,
    genesToMap = genesToMap, use = use, method = method
  )
  corMatch <- getTopMatch(corVar)
  colnames(corMatch) <- c("Class", "Correlation")

  dex <- apply(corVar, 1, function(x) return(diff(sort(-x)[1:2])))
  corMatch$DifferenceBetweenTopTwoCorrelations <- dex
  if(returnCor)
    corMatch <- cbind(corMatch,corVar)
  corMatch
}

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
    qprob = NULL,
    thresh = NULL,
    binarize = NULL,
    per_gene_correction = NULL,
    per_gene_layer = NULL,
    cell_type_col = NULL,
    area_col = NULL,
    x_col = NULL,
    y_col=NULL,
    cell_id= NULL
   
  )
  for (i in seq_along(args)) {
    if (args[i] %in% c("-s", "--spatial")) {
      arg_dict$spatial <- args[i + 1]
    } else if (args[i] %in% c("-d", "--dissociated")) {
      arg_dict$dissociated <- args[i + 1]
    }  else if (args[i] %in% c("-o", "--output")) {
      arg_dict$output <- args[i + 1]
    }  else if (args[i] %in% c("-p", "--hyperparamters")) {
      p <- args[i + 1]
    }  else if (args[i] %in% c("-g", "--group-paramters")) {
      g <- args[i + 1]
    }
  }

  # Check if -p parameter exists
  if (!is.null(p)) {
    # Parse JSON string
    json_string <- gsub("\\'", '"', p )  # Replace single quotes with double quotes
    json_string <- gsub(": True", ": true", json_string) 
    json_string <- gsub(": False", ": false", json_string) 
    params <- fromJSON(json_string)
    
    # Extract parameters
    if ("qprob" %in% names(params)) {
      arg_dict$qprob <- as.numeric(params$qprob)
    }
    if ("thresh" %in% names(params)) {
      arg_dict$thresh <- as.integer(params$thresh)
    }
    if ("binarize" %in% names(params)) {
      arg_dict$binarize <- as.logical(params$binarize)
    }
    if ("cell_type_col" %in% names(params)) {
      arg_dict$cell_type_col <- params$cell_type_col
    }
    if ("area_col" %in% names(params)) {
      arg_dict$area_col <- params$area_col
    }
    if ("x_col" %in% names(params)) {
      arg_dict$x_col <- params$x_col
    }
    if ("y_col" %in% names(params)) {
      arg_dict$y_col <- params$y_col
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
  hparams_defaults <- defaults$mfishtools 
  gparams_defaults <- defaults$annotation_params
  
  # Parse hyperparameters
  hyperparams <- list(
    qprob = if (!is.null(args$qprob)) args$qprob else hparams_defaults$qprob,
    thresh = if (!is.null(args$thresh)) args$thresh else hparams_defaults$thresh,
    binarize = if (!is.null(args$binarize)) args$binarize else hparams_defaults$binarize,
    cell_type_col = if (!is.null(args$cell_type_col)) args$cell_type_col else hparams_defaults$cell_type_col,
    area_col = if (!is.null(args$area_col)) args$area_col else hparams_defaults$area_col,
    x_col = if (!is.null(args$x_col)) args$x_col else hparams_defaults$x_col,
    y_col = if (!is.null(args$y_col)) args$y_col else hparams_defaults$y_col,
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
  spatial_data = t(as.matrix(adata$X)) # transpose matrix : genes as rows and cells as columns
  sc_data <- t(as.matrix(adata_sc$X)) 

 

  
  print('Preparing spatial data...')
  #selecting genes which appear in both spatial and sc data
  useGenes <- intersect(rownames(spatial_data),rownames(sc_data))
  spatial_data  <- spatial_data[useGenes,]
  # format metadata table
  selected_columns <- adata$obs[, c(hyperparams$area_col, hyperparams$x_col, hyperparams$y_col), drop = FALSE]

  # Define constant values for additional columns
  const_experiment <- "X"
  # Create new columns with constant values
  additional_columns <- data.frame(
    experiment = rep(const_experiment, nrow(adata$obs))
  )
  # Combine selected columns and additional columns to create a new data frame
  metadata <- cbind(selected_columns, additional_columns)
  names(metadata)[names(metadata) == hyperparams$x_col] <- "x"
  names(metadata)[names(metadata) == hyperparams$y_col] <- "y"


   # Cluster median calculation for single cell data
  print('Preparing Single cell data...')
  cl <- adata_sc$obs[[hyperparams$cell_type_col]] # Vector with cluster assignment of each cell (len = number of cells)


  # Initialize a list to store median expressions for each category
  medianExpr_list <- lapply(unique(cl), function(level) {
    columns_for_level <- sc_data[, cl == level, drop = FALSE]  # Subset based on matching categories
    rowMedians(columns_for_level)
  })
  # Combine median expressions into a single matrix
  medianExpr <- do.call("cbind", medianExpr_list)
  # Set column names to the unique categories
  colnames(medianExpr) <- unique(cl)
  rownames(medianExpr) <- rownames(sc_data)



  print("Mapping cell types...")
  #set weights
  weights <- NULL    # Integer weights.  SET TO NULL IF YOU DON'T KNOW WHAT YOU ARE DOING!
  #weights  <- round(rowSums(spatial_data)/min(rowSums(spatial_data)))  # Here is how to weight roughly by average expression level


  # Run cell type annotation
  log2p1   <- function(x) return(log2(x+1))  # log transform function
  fishMouse <- fishScaleAndMap(mapDat=spatial_data, refSummaryDat=medianExpr,
                            mappingFunction = cellToClusterMapping_byCor_, transform = log2p1, noiselevel = hyperparams$thresh, 
                            genesToMap = useGenes, metadata = metadata, qprob=hyperparams$qprob, binarize=hyperparams$binarize,
                            omitGenes = NULL,integerWeights=weights)


  annotation_df <- fishMouse$mappingResult
  annotation_df$cell_id <- adata$obs[[hyperparams$cell_id]]
  # Reorder columns in annotation_df with cell_id as the first column
  annotation_df <- annotation_df[, c("cell_id", setdiff(names(annotation_df), "cell_id"))]
  #rename columns
  colnames(annotation_df)[colnames(annotation_df) == "Class"] <- "celltype"
  colnames(annotation_df)[colnames(annotation_df) == "Correlation"] <- "score"
  
  annotation_df <- annotation_df[, -which(names(annotation_df) == "DifferenceBetweenTopTwoCorrelations")]

  # calculate probabilities from correlations
  calculate_probability <- function(row) {
    
    y <- row[5:length(row)]    # Extract values from the 5th column onwards
    scaledCorrelation <- pmax(y - (max(y) / 2), 0) ^ 2
    probability <- scaledCorrelation / sum(scaledCorrelation)
    return(probability)
  }


  probabilities <- t(apply(annotation_df[,5:ncol(annotation_df)], 1, calculate_probability))
  new_df <- cbind(annotation_df[, 1:4], probabilities)
  colnames(new_df) <- c(colnames(annotation_df)[1:4], colnames(probabilities))



  # Save annotation
  write.csv(new_df, file=args$output, row.names = FALSE)
}




# Entry point of the script
main <- function() {
  args <- parse_args()
  annotate_cells(args)
}

# Call main function
main()