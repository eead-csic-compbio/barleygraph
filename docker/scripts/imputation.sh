#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# --- 1. Argument Parsing ---
THREADS=8
READ_FILES=()

# Robustly parse flags and positional arguments in any order
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads) THREADS="$2"; shift 2 ;;
        *) READ_FILES+=("$1"); shift 1 ;;
    esac
done

if [ ${#READ_FILES[@]} -eq 0 ] || [ ${#READ_FILES[@]} -gt 2 ]; then
    echo "ERROR: Invalid number of read files."
    echo "Usage: $0 <read1.fq[.gz]> [read2.fq[.gz]] [-t <threads>]"
    exit 1
fi

READS_ARG=$(IFS=,; echo "${READ_FILES[*]}") # Joins array with commas
echo "Processing ${#READ_FILES[@]} read file(s) with ${THREADS} threads."
echo "Reads: ${READS_ARG}"

# --- 2. Project Detection & Variable Assignment ---
if [ -d "./Pan20/" ]; then
    PHG_PROJECT_DIR="./Pan20"
    PANGENOME_NAME="Pan20"
    
    echo "Pan20 found. Choose alignment method:"
    echo "1) gmap_geno"
    echo "2) mmap_pro"
    read -p "Enter choice (1 or 2): " VARIANT_CHOICE
    
    if [ "$VARIANT_CHOICE" == "1" ]; then VARIANT="gmap_geno";
    elif [ "$VARIANT_CHOICE" == "2" ]; then VARIANT="mmap_pro";
    else echo "ERROR: Invalid choice."; exit 1; fi

    HVCF_DIR="${PHG_PROJECT_DIR}/${VARIANT}/"
    OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/${VARIANT}"
    INDEX_PREFIX="${PANGENOME_NAME}_${VARIANT}"

elif [ -d "./Med13/" ]; then
    PHG_PROJECT_DIR="./Med13"
    PANGENOME_NAME="Med13"
    VARIANT=""
    HVCF_DIR="${PHG_PROJECT_DIR}/vcf_dbs/hvcf_files"
    OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/output"
    INDEX_PREFIX="${PANGENOME_NAME}"

elif [ -d "./Example_Ara/" ]; then
    PHG_PROJECT_DIR="./Example_Ara"
    PANGENOME_NAME="Example_Ara"
    VARIANT=""
    HVCF_DIR="${PHG_PROJECT_DIR}/vcf_dbs/hvcf_files"
    OUTPUT_BASE_DIR="${PHG_PROJECT_DIR}/output"
    INDEX_PREFIX="${PANGENOME_NAME}"

else
    echo "ERROR: Valid pangenome directory (Pan20, Med13, Example_Ara) not found."
    exit 1
fi

# Global Paths
IMPUTED_VCF_DIR="${OUTPUT_BASE_DIR}/imputed_vcf_files"
ROPEBWT_INDEX="${OUTPUT_BASE_DIR}/${INDEX_PREFIX}.fmd"
ASSEMBLIES_AGC="${PHG_PROJECT_DIR}/vcf_dbs/assemblies.agc"
SAMPLELIST="./samplelist.tsv"
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-/opt/conda/envs/phgv2.4/}"

mkdir -p "${IMPUTED_VCF_DIR}" "${PHG_PROJECT_DIR}/data"

echo "Pangenome: ${PANGENOME_NAME}${VARIANT:+ ($VARIANT)}"
echo "Index: ${ROPEBWT_INDEX}"
echo "--------------------------------------------------------"

# --- 3. Reference Genome Management ---
if [ ! -f "$SAMPLELIST" ]; then
    echo "ERROR: ${SAMPLELIST} not found."
    exit 1
fi

# Extract the first matching Reference genotype
REFERENCE_NAME=$(awk '$3=="Reference" {print $2; exit}' "$SAMPLELIST")
[[ -z "$REFERENCE_NAME" ]] && { echo "ERROR: No 'Reference' found in ${SAMPLELIST}."; exit 1; }

REFERENCE_GENOME="${PHG_PROJECT_DIR}/data/${REFERENCE_NAME}.fa"
echo "Detected Reference genome: ${REFERENCE_NAME}"

if [ ! -f "${ROPEBWT_INDEX}" ]; then
    echo "ERROR: Index missing at ${ROPEBWT_INDEX}. Run build_imputation_index.sh first."
    exit 1
fi

if [ ! -f "${REFERENCE_GENOME}" ]; then
    echo "Extracting ${REFERENCE_NAME} from AGC archive..."
    [[ ! -f "$ASSEMBLIES_AGC" ]] && { echo "ERROR: AGC archive missing at ${ASSEMBLIES_AGC}."; exit 1; }
    
    agc getset "$ASSEMBLIES_AGC" "$REFERENCE_NAME" > "${REFERENCE_GENOME}"
    
    if [ ! -s "${REFERENCE_GENOME}" ]; then
        echo "ERROR: Extraction failed. File is empty."
        rm -f "${REFERENCE_GENOME}" # Clean up empty file
        exit 1
    fi
fi

# --- 4. PHG Execution ---
export _JAVA_OPTIONS="-Xmx256g"

phg map-reads \
    --index "${ROPEBWT_INDEX}" \
    --read-files "${READS_ARG}" \
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

# --- 5. Cleanup ---
rm -f "${OUTPUT_BASE_DIR}/pathKeyFile.txt" "${REFERENCE_GENOME}"

echo "--------------------------------------------------------"
echo "PHG Imputation Completed. Output: ${IMPUTED_VCF_DIR}"
