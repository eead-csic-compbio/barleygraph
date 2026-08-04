#!/bin/bash

# Exit immediately if an unexpected error occurs
set -e

# ==============================================================================
# DATABASE CONFIGURATION
# To add a new database in the future, simply add a line in this format:
# "PANGENOME:SUBNAME:URL"
# ==============================================================================
DATABASES=(
    "Pan20:mmap_pro:http://example.com/path/to/Pan20-mmap_pro.tgz"
    "Pan20:gmap_geno:http://example.com/path/to/Pan20-gmap_geno.tgz"
    "Example_Ara:mmap_pro:http://example.com/path/to/Example_Ara-mmap_pro.tgz"
)

# Helper lookup maps
declare -A DB_URLS
declare -A DB_PANGENOMES

for entry in "${DATABASES[@]}"; do
    IFS=':' read -r pangenome subname url <<< "$entry"
    fullname="${pangenome}-${subname}"
    DB_URLS["$fullname"]="$url"
    DB_PANGENOMES["$fullname"]="$pangenome"
done

# Default database root
DATABASE_ROOT="${BARLEYGRAPH_DATABASE_PATH:-/barleygraph_databases}"
DATASET_NAME=""

# ------------------------------------------------------------------------------
# Intelligent Argument Parsing
# ------------------------------------------------------------------------------
if [ "$#" -gt 2 ]; then
  echo "Usage: $0 [dataset_name] [database_path]"
  exit 1
fi

for arg in "$@"; do
    # Remove trailing slashes from the argument just in case (e.g., "folder/" -> "folder")
    clean_arg="${arg%/}"
    
    if [[ -n "${DB_URLS[$arg]}" ]]; then
        # If the argument is a valid dataset name in our dictionary
        DATASET_NAME="$arg"
    else
        # If it's not a dataset name, assume it's the path
        DATABASE_ROOT="$clean_arg"
    fi
done

# ------------------------------------------------------------------------------
# Interactive Menu (Triggers if no valid dataset name was provided)
# ------------------------------------------------------------------------------
if [[ -z "$DATASET_NAME" ]]; then
  echo "Target directory set to: $DATABASE_ROOT"
  echo "Available databases:"
  echo ""

  # Extract unique pangenome names in order
  PANGENOMES=()
  for entry in "${DATABASES[@]}"; do
      IFS=':' read -r pg _ _ <<< "$entry"
      if [[ ! " ${PANGENOMES[*]} " =~ " ${pg} " ]]; then
          PANGENOMES+=("$pg")
      fi
  done

  # Print grouped menu
  declare -A MENU_MAP
  option_num=1

  for pg in "${PANGENOMES[@]}"; do
      echo "${pg}:"
      for entry in "${DATABASES[@]}"; do
          IFS=':' read -r p_name sub_name _ <<< "$entry"
          if [[ "$p_name" == "$pg" ]]; then
              echo "  ${option_num}- ${sub_name}"
              MENU_MAP["$option_num"]="${p_name}-${sub_name}"
              ((option_num++))
          fi
      done
      echo ""
  done

  echo -n "Enter the number of the database to setup: "
  read -r SELECTION

  if [[ -n "${MENU_MAP[$SELECTION]}" ]]; then
      DATASET_NAME="${MENU_MAP[$SELECTION]}"
  else
      echo "Invalid selection. Exiting."
      exit 1
  fi
fi

# ------------------------------------------------------------------------------
# Setup Execution
# ------------------------------------------------------------------------------
KIT_URL="${DB_URLS[$DATASET_NAME]}"
PANGENOME_NAME="${DB_PANGENOMES[$DATASET_NAME]}"

# Path definitions utilizing DATABASE_ROOT
PARENT_DIR="$DATABASE_ROOT/$PANGENOME_NAME"
TAR_FILE="$DATABASE_ROOT/${DATASET_NAME}.tgz"
EXTRACT_DIR="$PARENT_DIR/$DATASET_NAME"

echo "Using Database Root: $DATABASE_ROOT"

# Create the parent directory if it doesn't exist
if ! mkdir -p "$PARENT_DIR" 2>/dev/null; then
    echo "ERROR: Failed to create directory $PARENT_DIR."
    echo "Do you have write permissions? You may need to run this script with 'sudo'."
    exit 1
fi

# ------------------------------------------------------------------------------
# Download and Extract
# ------------------------------------------------------------------------------
# Create a hidden temporary directory for safe extraction (prevents naming collisions)
EXTRACT_DIR="$PARENT_DIR/.tmp_${DATASET_NAME}"

# Download archive if it doesn't exist in the root folder
if [ ! -f "$TAR_FILE" ]; then
  echo "Downloading ${DATASET_NAME}.tgz from $KIT_URL..."
  
  if ! wget -c -q --show-progress -O "$TAR_FILE" "$KIT_URL"; then
    echo "ERROR: Failed to download ${DATASET_NAME}.tgz."
    exit 1
  fi
fi

# Extract and merge files
if [ -f "$TAR_FILE" ]; then
  echo "Extracting $TAR_FILE..."
  
  mkdir -p "$EXTRACT_DIR"
  
  # Extract directly into the temporary folder
  if ! tar -xzf "$TAR_FILE" -C "$EXTRACT_DIR"; then
    echo "ERROR: Failed to extract ${DATASET_NAME}.tgz."
    exit 1
  fi

  # Merge contents into the parent folder (skipping existing/duplicate files)
  echo "Merging contents into shared $PANGENOME_NAME folder (skipping duplicate files)..."
  
  # cp -rn recursively merges folders and copies files without overwriting.
  cp -rn "$EXTRACT_DIR"/* "$PARENT_DIR"/ 2>/dev/null || true
  
  # Remove the temporary folder
  rm -rf "$EXTRACT_DIR"

  # Remove tarball after setup
  echo "Removing $TAR_FILE..."
  rm -f "$TAR_FILE"
fi

echo "Setup complete for $DATASET_NAME under $PANGENOME_NAME."
