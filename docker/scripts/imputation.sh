#!/bin/bash

# Exit immediately if an unexpected error occurs
set -e

# --- 0. Help Message ---
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <hvcf_folder> <read1.fq[.gz]> [read2.fq[.gz]] [-t <threads>] [--agc-file <path>]"
    echo ""
    echo "Description:"
    echo "  Runs PHG imputation by mapping reads against a RopeBWT index and finding paths."
    echo ""
    echo "Positional Arguments:"
    echo "  hvcf_folder       Path to the directory containing hVCF files and the .fmd index."
    echo "  read1.fq[.gz]     First FASTQ file (single-end or paired-end read 1)."
    echo "  read2.fq[.gz]     Second FASTQ file (optional, for paired-end)."
    echo ""
    echo "Options:"
    echo "  -t, --threads     Number of threads to use (default: 8)."
    echo "  --agc-file        Path to the AGC archive. If omitted, the script searches the hvcf_folder and its parent."
    echo "  -h, --help        Show this help message and exit."
    exit 0
fi

# --- 1. Argument Parsing ---
if [ "$#" -lt 2 ]; then
    echo "ERROR: Missing required arguments."
    echo "Run '$0 --help' for usage information."
    exit 1
fi

HVCF_DIR="${1%/}"
shift 1

THREADS=8
AGC_FILE=""
READ_FILES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads)
            THREADS="$2"
            if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
                echo "ERROR: Invalid value for threads."
                exit 1
            fi
            shift 2
            ;;
        --agc-file)
            AGC_FILE="$2"
            shift 2
            ;;
        -*)
            echo "ERROR: Unknown parameter $1"
            echo "Run '$0 --help' for usage information."
            exit 1
            ;;
        *)
            READ_FILES+=("$1")
            shift 1
            ;;
    esac
done

if [ ${#READ_FILES[@]} -eq 0 ] || [ ${#READ_FILES[@]} -gt 2 ]; then
    echo "ERROR: Invalid number of read files. Provide 1 (single-end) or 2 (paired-end)."
    exit 1
fi

# --- 2. Path Validation & Setup ---
if [ ! -d "$HVCF_DIR" ]; then
    echo "ERROR: hVCF directory '$HVCF_DIR' does not exist."
    exit 1
fi

# Determine the AGC file
if [ -z "$AGC_FILE" ]; then
    # Look in the hvcf directory using find to avoid 'set -e' crashing on ls failure
    AGC_FILE=$(find "${HVCF_DIR}" -maxdepth 1 -name "*.agc" | head -n 1)
    
    # If not found, check the parent directory
    if [ -z "$AGC_FILE" ]; then
        AGC_FILE=$(find "${HVCF_DIR}/.." -maxdepth 1 -name "*.agc" | head -n 1)
    fi
fi

if [ -z "$AGC_FILE" ] || [ ! -f "$AGC_FILE" ]; then
    echo "ERROR: Could not find an .agc file in '${HVCF_DIR}' or '${HVCF_DIR}/..'. Please provide one using --agc-file."
    exit 1
fi

# Automatically find the .fmd index built in the previous step
ROPEBWT_INDEX=$(find "${HVCF_DIR}" -maxdepth 1 -name "*.fmd" | head -n 1)

if [ -z "$ROPEBWT_INDEX" ] || [ ! -f "$ROPEBWT_INDEX" ]; then
    echo "ERROR: No .fmd index found in ${HVCF_DIR}. Run build_imputation_index first."
    exit 1
fi

READS_ARG=$(IFS=,; echo "${READ_FILES[*]}")
OUTPUT_BASE_DIR="${HVCF_DIR}/imputed_output"
IMPUTED_VCF_DIR="${OUTPUT_BASE_DIR}/imputed_vcf_files"
PATH_KEYFILE="${OUTPUT_BASE_DIR}/pathKeyFile.txt"

mkdir -p "${IMPUTED_VCF_DIR}"

# --- 3. AGC Reference Extraction ---
# Get the reference sample name dynamically from the archive
REF_NAME=$(agc listref "$AGC_FILE" | head -n 1)

if [ -z "$REF_NAME" ]; then
    echo "ERROR: Could not identify the reference sample inside $AGC_FILE."
    exit 1
fi

TEMP_FASTA="${OUTPUT_BASE_DIR}/${REF_NAME}.fa"

echo "Starting PHG Imputation"
echo "hVCF Dir:     ${HVCF_DIR}"
echo "Index File:   ${ROPEBWT_INDEX}"
echo "AGC Archive:  ${AGC_FILE}"
echo "Reference:    ${REF_NAME} (Extracting dynamically)"
echo "Reads:        ${READS_ARG}"
echo "Threads:      ${THREADS}"
echo "Output Dir:   ${IMPUTED_VCF_DIR}"
echo "--------------------------------------------------------"

echo "Extracting ${REF_NAME} to temporary FASTA..."
agc getset "$AGC_FILE" "$REF_NAME" > "$TEMP_FASTA"

if [ ! -s "$TEMP_FASTA" ]; then
    echo "ERROR: Extraction failed. The temporary FASTA file is empty."
    rm -f "$TEMP_FASTA"
    exit 1
fi

# --- 4. PHG Execution ---
export _JAVA_OPTIONS="-Xmx256g"
CONDA_ENV="/opt/miniconda3/envs/phgv2.4/"

echo "Running phg map-reads..."
phg map-reads \
    --index "${ROPEBWT_INDEX}" \
    --read-files "${READS_ARG}" \
    -o "${OUTPUT_BASE_DIR}" \
    --threads "${THREADS}" \
    --hvcf-dir "${HVCF_DIR}" \
    --conda-env-prefix "${CONDA_ENV}"

if [ ! -f "${PATH_KEYFILE}" ]; then
    echo "ERROR: phg map-reads failed to produce pathKeyFile.txt"
    exit 1
fi

echo "Running phg find-paths..."
phg find-paths \
    --path-keyfile "${PATH_KEYFILE}" \
    --hvcf-dir "${HVCF_DIR}" \
    --path-type haploid \
    --kmer-index "${ROPEBWT_INDEX}" \
    --threads "${THREADS}" \
    --reference-genome "${TEMP_FASTA}" \
    --output-dir "${IMPUTED_VCF_DIR}"

# --- 5. Cleanup ---
echo "Cleaning up temporary files..."
rm -f "${PATH_KEYFILE}" "${TEMP_FASTA}"

echo "--------------------------------------------------------"
echo "PHG Imputation Completed successfully."
echo "Imputed VCF files are located in: ${IMPUTED_VCF_DIR}"