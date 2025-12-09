#!/bin/bash

# --- Configuration Variables ---

# Parse command-line arguments
# Usage: ./imputation.sh ......

# Determine the base directory for the project
if [ -d "./Pan20/" ]; then
    PHG_PROJECT_DIR="./Pan20/"
elif [ -d "./Med13/" ]; then
    PHG_PROJECT_DIR="./Med13/"
else
    echo "ERROR: Neither ./Pan20/ nor ./Med13/ directories exist."
    exit 1
fi

# Store name of pangenome (folder name)
PANGENOME_NAME=$(basename "${PHG_PROJECT_DIR}")

# PHG Database Directories (Derived from project dir)
HVCF_DIR="${PHG_PROJECT_DIR}/vcf_dbs/hvcf_files"
OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/output"
IMPUTED_VCF_DIR="${OUTPUT_BASE_DIR}/imputed_vcf_files"

# Ensure output directories exist
mkdir -p "${OUTPUT_BASE_DIR}"   # Actually redundant, but safe
mkdir -p "${IMPUTED_VCF_DIR}"

# Pre-built Index and Reference (You should update these paths)
ROPEBWT_INDEX="${OUTPUT_BASE_DIR}/${PANGENOME_NAME}.fmd"
REFERENCE_GENOME="${PHG_PROJECT_DIR}/data/MorexV3.fa"

# Ensure required files exist

if [ ! -f "${ROPEBWT_INDEX}" ]; then
    echo "ERROR: RopeBWT index not found at ${ROPEBWT_INDEX}. Please build the index first."
    echo "Run build_imputation_index.sh script."
    exit 1
fi

if [ ! -f "${REFERENCE_GENOME}" ]; then
    echo "ERROR: Reference genome not found at ${REFERENCE_GENOME}."
    echo "Extracting it and saving to ${REFERENCE_GENOME}..."
    agc getset "${PHG_PROJECT_DIR}/vcf_dbs/assemblies.agc" MorexV3 > "${REFERENCE_GENOME}"
fi

# --- Main Script Execution ---

# 1. Input Validation and Variable Assignment
if [ $# -lt 2 ] || [ $# -gt 4 ]; then
    echo "Usage: $0 <read1.fq[.gz]> [read2.fq[.gz]] -t <threads>"
    echo "Example (single-end): $0 reads/sample_r1.fq.gz -t 8"
    echo "Example (paired-end): $0 reads/sample_r1.fq.gz reads/sample_r2.fq.gz -t 8"
    exit 1
fi

READ_FILE_1="$1" # This is always present

if [ $# -eq 4 ]; then   # If 4 arguments mean 2 fastq and threads specified
    READ_FILE_2="$2"    # So this is the second read file
    if [ "$3" != "-t" ]; then
        echo "ERROR: Third argument must be '-t'."
        exit 1
    else
        THREADS="$4"    # Assign the number of threads
    fi

    echo "Processing paired-end reads with ${THREADS} threads."
    echo "Read 1: ${READ_FILE_1}"
    echo "Read 2: ${READ_FILE_2}"

elif [ $# -eq 3 ]; then   # If 3 arguments, it must be 1 read file and threads specified
    if [ "$2" != "-t" ]; then
        echo "ERROR: Second argument must be '-t'."
        exit 1
    else
        THREADS="$3"    # Assign the number of threads
    fi
    echo "Processing single-end reads with ${THREADS} threads."
    echo "Read 1: ${READ_FILE_1}"
    READ_FILE_2=""      # No second read file

elif [ $# -eq 2 ]; then   # If only 2 arguments, it must be 2 read files without threads
    READ_FILE_2="$2"    # So this is the second read file
    THREADS=8            # Default threads
    echo "Processing paired-end reads with default 8 threads."
    echo "Read 1: ${READ_FILE_1}"
    echo "Read 2: ${READ_FILE_2}"

elif [ $# -eq 1 ]; then   # If only 1 argument, it must be 1 read file without threads
    THREADS=8            # Default threads
    echo "Processing single-end reads with default 8 threads."
    echo "Read 1: ${READ_FILE_1}"
    READ_FILE_2=""      # No second read file

else
    echo "ERROR: Invalid number of arguments."
    echo "Usage: $0 <read1.fq[.gz]> [read2.fq[.gz]] [-t <threads>]"
    exit 1
fi

echo "Starting PHG pipeline"
echo "--------------------------------------------------------"



# Increase ram to avoid Java heap space issues
export _JAVA_OPTIONS="-Xmx256g"

# Construct read-files argument (comma-separated for paired-end)
if [ -n "${READ_FILE_2}" ]; then
    READ_FILES="${READ_FILE_1},${READ_FILE_2}"
else
    READ_FILES="${READ_FILE_1}"
fi

# Allow overriding conda prefix via env var; fallback to original default
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-/opt/conda/envs/phgv2.4/}"

phg map-reads --index "${ROPEBWT_INDEX}" \
    --read-files "${READ_FILES}" \
    -o "${OUTPUT_BASE_DIR}" \
    --threads "${THREADS}" \
    --hvcf-dir "${HVCF_DIR}" \
    --conda-env-prefix "${CONDA_ENV_PREFIX}"

phg find-paths \
    --path-keyfile "${OUTPUT_BASE_DIR}/pathKeyFile.txt" \
    --hvcf-dir "${HVCF_DIR}" \
    --path-type haploid \
    --kmer-index "${ROPEBWT_INDEX}" \
    --threads "${THREADS}" \
    --reference-genome "${REFERENCE_GENOME}" \
    --output-dir "${IMPUTED_VCF_DIR}"

echo "PHG Imputation Completed. Imputed VCF files are located in: ${IMPUTED_VCF_DIR}"
echo "--------------------------------------------------------"

# Remove intermediate files to save space
rm -rf "${OUTPUT_BASE_DIR}/pathKeyFile.txt"
rm -f "${REFERENCE_GENOME}"
