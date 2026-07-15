#!/bin/bash

# --- Configuration Variables ---

# Parse command-line arguments
# Usage: ./imputation.sh ......

# Determine the base directory for the project
if [ -d "./Pan20/" ]; then
    PHG_PROJECT_DIR="./Pan20/"
    PANGENOME_NAME="Pan20"
    
    # For Pan20, ask which variant to use
    echo "Pan20 found. Choose alignment method:"
    echo "1) gmap-geno"
    echo "2) mmap-pro"
    read -p "Enter choice (1 or 2): " VARIANT_CHOICE
    
    if [ "$VARIANT_CHOICE" == "1" ]; then
        VARIANT="gmap-geno"
    elif [ "$VARIANT_CHOICE" == "2" ]; then
        VARIANT="mmap-pro"
    else
        echo "ERROR: Invalid choice. Please enter 1 or 2."
        exit 1
    fi
    
elif [ -d "./Med13/" ]; then
    PHG_PROJECT_DIR="./Med13/"
    PANGENOME_NAME="Med13"
    VARIANT=""

elif [ -d "./Example_Ara/" ]; then
    PHG_PROJECT_DIR="./Example_Ara/"
    PANGENOME_NAME="Example_Ara"
    VARIANT=""

else
    echo "ERROR: Neither ./Pan20/ nor ./Med13/ nor ./Example_Ara/ directories exist."
    exit 1
fi

# PHG Database Directories (Derived from project dir)
HVCF_DIR="${PHG_PROJECT_DIR}/vcf_dbs/hvcf_files"

# For Pan20 variants, use variant-specific output directories
if [ -n "$VARIANT" ]; then
    OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/${VARIANT}"
    INDEX_PREFIX="${PANGENOME_NAME}_${VARIANT}"
else
    OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/output"
    INDEX_PREFIX="${PANGENOME_NAME}"
fi

IMPUTED_VCF_DIR="${OUTPUT_BASE_DIR}/imputed_vcf_files"

# Ensure output directories exist
mkdir -p "${OUTPUT_BASE_DIR}"   # Actually redundant, but safe
mkdir -p "${IMPUTED_VCF_DIR}"

# Extract the Reference genotype name dynamically from samplelist.tsv
SAMPLELIST="./samplelist.tsv"
if [ ! -f "$SAMPLELIST" ]; then
    echo "ERROR: ${SAMPLELIST} not found. Cannot determine Reference genome."
    exit 1
fi

REFERENCE_NAME=$(awk '$3=="Reference" {print $2}' "$SAMPLELIST" | head -n 1)

if [ -z "$REFERENCE_NAME" ]; then
    echo "ERROR: Could not find a genotype marked as 'Reference' in ${SAMPLELIST}."
    exit 1
fi

echo "Detected Reference genome: ${REFERENCE_NAME}"

# Pre-built Index and Reference
ROPEBWT_INDEX="${OUTPUT_BASE_DIR}/${INDEX_PREFIX}.fmd"
REFERENCE_GENOME="${PHG_PROJECT_DIR}/data/${REFERENCE_NAME}.fa"

# Ensure required files exist
if [ ! -f "${ROPEBWT_INDEX}" ]; then
    echo "ERROR: RopeBWT index not found at ${ROPEBWT_INDEX}. Please build the index first."
    echo "Run build_imputation_index.sh script."
    exit 1
fi

if [ ! -f "${REFERENCE_GENOME}" ]; then
    echo "ERROR: Reference genome not found at ${REFERENCE_GENOME}."
    echo "Extracting ${REFERENCE_NAME} from AGC archive and saving to ${REFERENCE_GENOME}..."
    agc getset "${PHG_PROJECT_DIR}/vcf_dbs/assemblies.agc" "${REFERENCE_NAME}" > "${REFERENCE_GENOME}"
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
echo "Pangenome: ${PANGENOME_NAME}${VARIANT:+ ($VARIANT)}"
echo "Index: ${ROPEBWT_INDEX}"
echo "Imputed VCF output: ${IMPUTED_VCF_DIR}"
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
echo "Index used: ${ROPEBWT_INDEX}"
echo "--------------------------------------------------------"

# Remove intermediate files to save space
rm -rf "${OUTPUT_BASE_DIR}/pathKeyFile.txt"
rm -f "${REFERENCE_GENOME}"
