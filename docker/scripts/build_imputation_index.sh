#!/bin/bash

# --- Configuration Variables ---

# Determine the base directory for the project
if [ -d "./Pan20/" ]; then
    PHG_PROJECT_DIR="./Pan20/"
elif [ -d "./Med13/" ]; then
    PHG_PROJECT_DIR="./Med13/"
else
    echo "ERROR: Neither ./Pan/ nor ./Med13/ directories exist."
    exit 1
fi
# Store name of pangenome (folder name)
PANGENOME_NAME=$(basename "${PHG_PROJECT_DIR}")

# PHG Database Directories (Derived from project dir)
DB_PATH="${PHG_PROJECT_DIR}/vcf_dbs"
OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/output"
INDEX_OUTPUT_DIR="${OUTPUT_BASE_DIR}"

# Define the full path for the output index file
ROPEBWT_INDEX="${INDEX_OUTPUT_DIR}/${PANGENOME_NAME}.fmd"

# Default settings
DEFAULT_THREADS=8
THREADS=${DEFAULT_THREADS}

# Ensure output directory exists
mkdir -p "${INDEX_OUTPUT_DIR}"

# --- Main Script Execution ---

# Check if it is already built and found in gmap_db or at output location
if [ -f "${DB_PATH}/gmap_db/$(basename ${ROPEBWT_INDEX})" ] || [ -f "${ROPEBWT_INDEX}" ]; then
    echo "PHG RopeBWT Index already exists at:"
    if [ -f "${DB_PATH}/gmap_db/$(basename ${ROPEBWT_INDEX})" ]; then
        echo "${DB_PATH}/gmap_db/$(basename ${ROPEBWT_INDEX})"
    else
        echo "${ROPEBWT_INDEX}"
    fi

    # generate the symlink if not already present
    if [ ! -f "${ROPEBWT_INDEX}" ]; then
        ln -s ${DB_PATH}/gmap_db/$(basename ${ROPEBWT_INDEX}) ${ROPEBWT_INDEX}
        echo "Created symlink at: ${ROPEBWT_INDEX}"
    fi

    echo "Skipping index build."
    exit 0
fi

# 1. Input Validation and Argument Assignment
if [ $# -eq 0 ]; then
    # No arguments provided, use default threads
    echo "No threads specified. Using default THREADS=${DEFAULT_THREADS}."

elif [ $# -eq 2 ]; then
    # Two arguments: -t <threads>
    if [ "$1" == "-t" ]; then
        THREADS="$2"
        # Basic check to ensure threads is a number
        if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
            echo "ERROR: Invalid value for threads. Please provide a positive integer."
            echo "Usage: $0 [-t <threads>]"
            exit 1
        fi
    else
        echo "ERROR: Invalid argument format."
        echo "Usage: $0 [-t <threads>]"
        exit 1
    fi
else
    echo "ERROR: Invalid number of arguments."
    echo "Usage: $0 [-t <threads>]"
    exit 1
fi


echo "Starting PHG RopeBWT Index Build for Pangenome: ${PANGENOME_NAME}"
echo "Threads: ${THREADS}"
echo "DB Path: ${DB_PATH}"
echo "Output Index: ${ROPEBWT_INDEX}"
echo "--------------------------------------------------------"

# Increase ram to avoid Java heap space issues
export _JAVA_OPTIONS="-Xmx256g"

# 2. Execute phg rope-bwt-index
phg rope-bwt-index \
    --db-path "${DB_PATH}" \
    --hvcf-dir "${DB_PATH}/hvcf_files/" \
    --output-dir "${INDEX_OUTPUT_DIR}" \
    --index-file-prefix "${PANGENOME_NAME}" \
    --threads "${THREADS}" \
    --delete-fmr-index \
    --conda-env-prefix /opt/conda/envs/phgv2.4/

# 3. Check for successful completion
if [ $? -ne 0 ]; then
    echo "ERROR: phg rope-bwt-index failed. Check the logs above."
    exit 1
fi

# 4. Final verification
if [ -f "${ROPEBWT_INDEX}" ]; then
    echo "--------------------------------------------------------"
    echo "PHG Index built successfully."
    # Move it to gmap_db folder to keep it for future use
    mv ${ROPEBWT_INDEX} ${DB_PATH}/gmap_db/
    echo "Moved index to: ${DB_PATH}/gmap_db/"
    # but keep a direct link to use it
    ln -s ${DB_PATH}/gmap_db/$(basename ${ROPEBWT_INDEX}) ${ROPEBWT_INDEX}
    echo "Index link created at:"
    echo "${ROPEBWT_INDEX}"
    echo "Removed intermediate files."
    echo "--------------------------------------------------------"
else
    echo "ERROR: Index file not found after running the command. Build failed."
    exit 1
fi
