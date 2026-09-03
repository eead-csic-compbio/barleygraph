#!/bin/bash

# Exit immediately if an unexpected error occurs
set -e

# --- 0. Help Message ---
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <hvcf_folder> [-t <threads>] [--db-path <path>]"
    echo ""
    echo "Description:"
    echo "  Builds a PHG RopeBWT index (.fmd) from a directory of hVCF files."
    echo "  The index will be saved directly into the provided hVCF folder."
    echo ""
    echo "Positional Arguments:"
    echo "  hvcf_folder       Path to the directory containing hVCF files."
    echo ""
    echo "Options:"
    echo "  -t, --threads     Number of threads to use (default: 8)."
    echo "  --db-path         Path to the folder containing the .agc file."
    echo "                    If omitted, the script searches the hvcf_folder and its parent."
    echo "  -h, --help        Show this help message and exit."
    exit 0
fi

# 1. Require the hVCF folder as the first positional argument
if [ -z "$1" ] || [[ "$1" == "-"* ]]; then
    echo "ERROR: You must provide the path to the folder containing hVCF files."
    echo "Run '$0 --help' for usage information."
    exit 1
fi

HVCF_DIR="${1%/}"

if [ ! -d "$HVCF_DIR" ]; then
    echo "ERROR: hVCF directory '$HVCF_DIR' does not exist."
    exit 1
fi

shift # Move past the folder argument

# Default settings
THREADS=8
DB_PATH=""

# Parse remaining optional arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads)
            THREADS="$2"
            if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
                echo "ERROR: Invalid value for threads. Please provide a positive integer."
                exit 1
            fi
            shift 2
            ;;
        --db-path)
            DB_PATH="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown parameter: $1"
            echo "Run '$0 --help' for usage information."
            exit 1
            ;;
    esac
done

# --- DB Path Auto-Discovery ---
if [ -z "$DB_PATH" ]; then
    # Check if an .agc file exists in the hVCF directory
    if [ -n "$(find "${HVCF_DIR}" -maxdepth 1 -name "*.agc" | head -n 1)" ]; then
        DB_PATH="${HVCF_DIR}"
    # If not found, check the parent directory
    elif [ -n "$(find "${HVCF_DIR}/.." -maxdepth 1 -name "*.agc" | head -n 1)" ]; then
        DB_PATH="${HVCF_DIR}/.."
    fi
fi

if [ -z "$DB_PATH" ] || [ ! -d "$DB_PATH" ]; then
    echo "ERROR: Could not find a folder containing an .agc file automatically."
    echo "Please manually provide the path using: --db-path <folder_with_agc>"
    exit 1
fi

# Convert DB_PATH to absolute path to avoid issues if relative path changes
DB_PATH=$(realpath "$DB_PATH")

# --- Set Paths ---
PANGENOME_NAME=$(basename "$HVCF_DIR")

# Set output folder to be EXACTLY the hVCF folder
INDEX_OUTPUT_DIR="${HVCF_DIR}"
INDEX_PREFIX="${PANGENOME_NAME}"
ROPEBWT_INDEX="${INDEX_OUTPUT_DIR}/${INDEX_PREFIX}.fmd"

# --- Main Script Execution ---

if [ -f "${ROPEBWT_INDEX}" ]; then
    echo "PHG RopeBWT Index already exists at: ${ROPEBWT_INDEX}"
    echo "Skipping index build."
    exit 0
fi

echo "Starting PHG RopeBWT Index Build"
echo "hVCF Dir:     ${HVCF_DIR}"
echo "DB Path:      ${DB_PATH}"
echo "Output Index: ${ROPEBWT_INDEX}"
echo "Threads:      ${THREADS}"
echo "--------------------------------------------------------"

# Increase ram to avoid Java heap space issues
export _JAVA_OPTIONS="-Xmx256g"

# Execute phg rope-bwt-index
phg rope-bwt-index \
    --db-path "${DB_PATH}" \
    --hvcf-dir "${HVCF_DIR}" \
    --output-dir "${INDEX_OUTPUT_DIR}" \
    --index-file-prefix "${INDEX_PREFIX}" \
    --threads "${THREADS}" \
    --delete-fmr-index \
    --conda-env-prefix /opt/miniconda3/envs/phgv2.4/

if [ $? -ne 0 ]; then
    echo "ERROR: phg rope-bwt-index failed. Check the logs above."
    exit 1
fi

if [ -f "${ROPEBWT_INDEX}" ]; then
    echo "--------------------------------------------------------"
    echo "PHG Index built successfully at:"
    echo "${ROPEBWT_INDEX}"
    echo "--------------------------------------------------------"
# --- Cleanup ---
    # Remove the intermediate pangenome.fa file
    rm -f "${INDEX_OUTPUT_DIR}/pangenome.fa" "pangenome.fa"
    
else
    echo "ERROR: Index file not found after running the command. Build failed."
    exit 1
fi