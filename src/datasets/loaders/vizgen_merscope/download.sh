#!/bin/bash

######################################################################
# Requirements: 
# - Install google-cloud-sdk, e.g. mamba install -c conda-forge google-cloud-sdk
# - Authenticate with your Google account: gcloud auth login --no-browser
#
# Usage: merfish_download.sh <BUCKET_NAME> <DATASET> <OUT_DIR> <DRY_RUN=true/false>
# - DATASET refers to the naming in the Vizgen download page
#
# Vizgen cancer panel datasets download page: https://console.cloud.google.com/storage/browser/vz-ffpe-showcase
# The script should also work for other Vizgen style data buckets.

BUCKET_NAME=$1
DATASET=$2
OUT_DIR=$3
DRY_RUN=$4

######################################################################

# Suppress BrokenPipeError messages from Python
export PYTHONWARNINGS="ignore:BrokenPipeError"

LIST_LIMIT=10  # Maximum number of files to list in dry-run mode

# List of files or directories to download
FILES=(
    "cell_by_gene.csv"
    "cell_metadata.csv"
    "detected_transcripts.csv"
    "cell_boundaries/"  # This is a directory
    "cell_boundaries.parquet"
    "images/mosaic_DAPI_z3.tif" # If you want all images, change this to "images/"
)

# Check if gsutil is installed
if ! command -v gsutil &> /dev/null; then
    echo "ERROR: gsutil is not installed. Please install the Google Cloud SDK."
fi

# Check if user is authenticated
if ! gcloud auth list 2>&1 | grep -q ACTIVE; then
    echo "ERROR: No authenticated Google Cloud account found. Please authenticate using 'gcloud auth login'."
fi


# Loop over each file or directory to download
for FILE in "${FILES[@]}"; do
    # Get the directory from the file path (remove the file name)
    DIR=$(dirname "$FILE")

    # Create the full target directory path
    TARGET_DIR="$OUT_DIR/$DIR"

    # If DIR is '.', set TARGET_DIR to OUT_DIR
    if [ "$DIR" = "." ]; then
        TARGET_DIR="$OUT_DIR"
    fi

    # Create the directory if it doesn't exist
    mkdir -p "$TARGET_DIR"
    
    # Check if the file or directory exists in the bucket
    if gsutil ls -d "gs://$BUCKET_NAME/$DATASET/$FILE" > /dev/null 2>&1; then
        if [ "$DRY_RUN" = true ]; then
            # Dry-run: list files instead of downloading
            echo -e "\nSimulating download for gs://$BUCKET_NAME/$DATASET/$FILE"
            echo -e "Download to $TARGET_DIR"

            if [[ "$FILE" == */ ]]; then
                echo "Listing first $LIST_LIMIT files from gs://$BUCKET_NAME/$DATASET/$FILE"
                gsutil ls "gs://$BUCKET_NAME/$DATASET/$FILE" 2>/dev/null | head -n $LIST_LIMIT
            else
                gsutil ls "gs://$BUCKET_NAME/$DATASET/$FILE"
            fi
        else
            # If it's a directory, download recursively
            if [[ "$FILE" == */ ]]; then
                echo "Downloading directory gs://$BUCKET_NAME/$DATASET/$FILE to $TARGET_DIR/"
                gsutil -m cp -r "gs://$BUCKET_NAME/$DATASET/$FILE" "$TARGET_DIR/"
            else
                # Download the file, maintaining directory structure under OUT_DIR
                echo "Downloading file gs://$BUCKET_NAME/$DATASET/$FILE to $TARGET_DIR/"
                gsutil cp "gs://$BUCKET_NAME/$DATASET/$FILE" "$TARGET_DIR/"
            fi
        fi
    else
        # Skip download if file or directory does not exist
        echo "Skipping $FILE: does not exist in gs://$BUCKET_NAME/$DATASET/"
    fi
done
