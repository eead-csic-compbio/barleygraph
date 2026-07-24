#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build_PHG_database.sh --config <path/to/database.config> [--dry-run]

The config file must define at least:
  DATABASE_NAME=<name>
  REFERENCE_FASTA=<path-or-url>
  REFERENCE_GFF=<path-or-url>

Add one assembly source per line, for example:
  ASSEMBLYS=/path/to/genome1.fa
  ASSEMBLYS=/path/to/genome2.fa
  ASSEMBLYS=https://example.org/genome3.fa

Optional keys:
  REFERENCE_NAME=<name>          # defaults to basename of REFERENCE_FASTA
  THREADS=<int>                  # default 32
  CONDA_ENV=<name>               # default phgv2-conda
  GMAP_BUILD_CMD=<path>          # default /usr/local/bin/gmap_build
EOF
}

CONFIG_FILE=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)
            [[ $# -ge 2 ]] || { echo "Error: missing value for $1" >&2; exit 1; }
            CONFIG_FILE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$CONFIG_FILE" ]]; then
    echo "Error: --config is required." >&2
    usage >&2
    exit 1
fi

check_requirements() {
    local missing=()

    if ! command -v tabix >/dev/null 2>&1; then
        missing+=("tabix (install with conda install -c bioconda tabix, apt-get install tabix, or brew install tabix)")
    fi

    if ! python3 - <<'PY' >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec('tqdm') else 1)
PY
    then
        missing+=("tqdm (install with pip install tqdm)")
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "Missing required dependencies:" >&2
        for dep in "${missing[@]}"; do
            echo " - $dep" >&2
        done
        echo "Please install the dependencies above and rerun the script." >&2
        exit 1
    fi
}

check_requirements

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config file not found: $CONFIG_FILE" >&2
    exit 1
fi

if [[ "$CONFIG_FILE" != *.config ]]; then
    echo "Error: config file must have a .config extension: $CONFIG_FILE" >&2
    exit 1
fi

CONFIG_ABS="$CONFIG_FILE"
if [[ "$CONFIG_FILE" != /* ]]; then
    CONFIG_ABS="$PWD/$CONFIG_FILE"
fi
CONFIG_DIR="$(cd "$(dirname "$CONFIG_ABS")" && pwd)"
CONFIG_FILE_ABS="$CONFIG_DIR/$(basename "$CONFIG_ABS")"

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

strip_quotes() {
    local value="$1"
    value="${value#\"}"
    value="${value%\"}"
    value="${value#\'}"
    value="${value%\'}"
    printf '%s' "$value"
}

split_config_values() {
    local value="$1"
    local token
    while IFS= read -r token; do
        token="$(trim "$token")"
        if [[ -n "$token" ]]; then
            printf '%s\n' "$token"
        fi
    done < <(printf '%s\n' "$value" | tr ',;' '\n')
}

append_assembly_sources() {
    local value="$1"
    local item
    while IFS= read -r item; do
        item="$(trim "$item")"
        [[ -n "$item" ]] && ASSEMBLY_SOURCES+=("$item")
    done < <(split_config_values "$value")
}

is_url() {
    local value="$1"
    [[ "$value" =~ ^https?:// ]] || [[ "$value" =~ ^ftp:// ]]
}

resolve_path() {
    local value="$1"
    if [[ "$value" =~ ^/ ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$CONFIG_DIR/$value"
    fi
}

basename_without_ext() {
    local path="$1"
    local name
    name="$(basename "$path")"
    name="${name%.gz}"
    name="${name%.fasta}"
    name="${name%.fa}"
    name="${name%.fna}"
    printf '%s\n' "$name"
}

resolve_gmap_build_cmd() {
    local cmd="$1"
    local base
    base="$(basename "$cmd")"

    if [[ "$base" == "gmap" ]]; then
        if [[ "$cmd" == */* ]]; then
            local dir
            dir="$(dirname "$cmd")"
            if [[ -x "$dir/gmap_build" ]]; then
                printf '%s\n' "$dir/gmap_build"
                return 0
            fi
        fi

        if command -v gmap_build >/dev/null 2>&1; then
            command -v gmap_build
            return 0
        fi
    fi

    printf '%s\n' "$cmd"
}

run_cmd() {
    local -a cmd=("$@")
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '[DRY-RUN] '
        printf '%q ' "${cmd[@]}"
        echo
    else
        "${cmd[@]}"
    fi
}

require_command() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: required command not found: $cmd" >&2
        exit 1
    fi
}

require_python_module() {
    local module="$1"
    local install_hint="$2"

    if ! python3 - <<'PY' "$module" >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec(sys.argv[1]) else 1)
PY
    then
        echo "Error: required Python package not found: $module" >&2
        echo "Please install it with: $install_hint" >&2
        exit 1
    fi
}

check_requirements() {
    local missing=()

    if ! command -v tabix >/dev/null 2>&1; then
        missing+=("tabix (install with conda install -c bioconda tabix, apt-get install tabix, or brew install tabix)")
    fi

    if ! python3 - <<'PY' >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec('tqdm') else 1)
PY
    then
        missing+=("tqdm (install with pip install tqdm)")
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "Missing required dependencies:" >&2
        for dep in "${missing[@]}"; do
            echo " - $dep" >&2
        done
        echo "Please install the dependencies above and rerun the script." >&2
        exit 1
    fi
}

index_vcf_files() {
    local dir
    local file

    if ! command -v tabix >/dev/null 2>&1; then
        echo "Error: tabix is required to index VCF files" >&2
        exit 1
    fi

    for dir in "$@"; do
        [[ -d "$dir" ]] || continue
        while IFS= read -r file; do
            [[ -n "$file" ]] || continue
            if [[ ! -f "$file.tbi" && ! -f "$file.csi" ]]; then
                if [[ "$DRY_RUN" -eq 1 ]]; then
                    echo "[DRY-RUN] would index $file"
                else
                    echo "Indexing $file"
                    tabix -f -p vcf "$file"
                fi
            fi
        done < <(find "$dir" -maxdepth 1 -type f \( -name '*.vcf.gz' -o -name '*.vcf' \) | sort)
    done
}

activate_conda() {
    if [[ -n "${CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
        if [[ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1090
            source "$(conda info --base)/etc/profile.d/conda.sh"
        fi
        conda activate "$CONDA_ENV" >/dev/null 2>&1 || true
    fi
}

prepare_source_file() {
    local source="$1"
    local sample_name="$2"
    local target_file="$3"
    local tmp_download=""

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] would prepare source $source -> $target_file"
        return 0
    fi

    if [[ -f "$target_file" ]]; then
        echo "Using existing prepared file: $target_file"
        return 0
    fi

    if is_url "$source"; then
        tmp_download="$(mktemp)"
        echo "Downloading $source"
        if command -v curl >/dev/null 2>&1; then
            curl -L --fail -o "$tmp_download" "$source"
        elif command -v wget >/dev/null 2>&1; then
            wget -O "$tmp_download" "$source"
        else
            echo "Error: neither curl nor wget is available for URL downloads" >&2
            exit 1
        fi
    else
        local source_abs
        source_abs="$(resolve_path "$source")"
        if [[ ! -f "$source_abs" ]]; then
            echo "Error: source file not found: $source_abs" >&2
            exit 1
        fi
        tmp_download="$source_abs"
    fi

    if [[ "$tmp_download" == *.gz ]] || gzip -t "$tmp_download" >/dev/null 2>&1; then
        gzip -dc "$tmp_download" > "$target_file"
    else
        cp "$tmp_download" "$target_file"
    fi

    python3 - "$target_file" "$sample_name" <<'PY'
import sys
from pathlib import Path
import re

path = Path(sys.argv[1])
sample_name = sys.argv[2]
lines = path.read_text().splitlines()
new_lines = []
for line in lines:
    if line.startswith('>'):
        if re.match(r'^>ENA\|', line):
            m = re.search(r'chromosome:\s*(\S+)', line)
            if m:
                new_lines.append(f">chr{m.group(1)} sampleName={sample_name}")
            else:
                new_lines.append(f">{line[1:]} sampleName={sample_name}")
        elif re.match(r'^>chr', line):
            new_lines.append(f">{line[1:]} sampleName={sample_name}")
        else:
            new_lines.append(f">{line[1:]} sampleName={sample_name}")
    else:
        new_lines.append(line)
path.write_text("\n".join(new_lines) + "\n")
PY
}

# Parse config
DATABASE_NAME=""
REFERENCE_FASTA=""
REFERENCE_GFF=""
REFERENCE_NAME=""
THREADS="32"
CONDA_ENV="phgv2-conda"
GMAP_BUILD_CMD="/usr/local/bin/gmap_build"
ASSEMBLY_SOURCES=()

while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
    line="${raw_line%%$'\r'}"
    line="${line%%#*}"
    line="$(trim "$line")"
    [[ -z "$line" ]] && continue

    if [[ "$line" =~ ^([A-Za-z0-9_]+)[[:space:]]*=(.*)$ ]]; then
        key="${BASH_REMATCH[1]}"
        value="$(strip_quotes "$(trim "${BASH_REMATCH[2]}")")"
        case "$key" in
            DATABASE_NAME)
                DATABASE_NAME="$value"
                ;;
            REFERENCE_FASTA)
                REFERENCE_FASTA="$value"
                ;;
            REFERENCE_GFF)
                REFERENCE_GFF="$value"
                ;;
            REFERENCE_NAME)
                REFERENCE_NAME="$value"
                ;;
            THREADS)
                THREADS="$value"
                ;;
            CONDA_ENV)
                CONDA_ENV="$value"
                ;;
            GMAP_BUILD_CMD)
                GMAP_BUILD_CMD="$value"
                ;;
            ASSEMBLYS|ASSEMBLY_FASTA|ASSEMBLY_SOURCE|ASSEMBLY_FASTAS)
                append_assembly_sources "$value"
                ;;
            ASSEMBLY_FASTA_*)
                ASSEMBLY_SOURCES+=("$value")
                ;;
        esac
    fi
done < "$CONFIG_FILE_ABS"

[[ -n "$DATABASE_NAME" ]] || { echo "Error: DATABASE_NAME is required in the config file" >&2; exit 1; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: REFERENCE_FASTA is required in the config file" >&2; exit 1; }
[[ -n "$REFERENCE_GFF" ]] || { echo "Error: REFERENCE_GFF is required in the config file" >&2; exit 1; }
[[ ${#ASSEMBLY_SOURCES[@]} -gt 0 ]] || { echo "Error: at least one assembly source must be provided via ASSEMBLYS or ASSEMBLY_FASTA_N" >&2; exit 1; }

ROOT_DIR="$CONFIG_DIR/$DATABASE_NAME"
GMAP_BUILD_CMD="$(resolve_gmap_build_cmd "$GMAP_BUILD_CMD")"

if [[ -z "$REFERENCE_NAME" ]]; then
    REFERENCE_NAME="$(basename_without_ext "$REFERENCE_FASTA")"
fi
DATA_DIR="$ROOT_DIR/data"
OUTPUT_DIR="$ROOT_DIR/output"
VCF_DB_DIR="$ROOT_DIR/vcf_dbs"
GMAP_DB_DIR="$ROOT_DIR/gmap_db"

mkdir -p "$ROOT_DIR" "$DATA_DIR" "$OUTPUT_DIR" "$VCF_DB_DIR" "$OUTPUT_DIR/alignment_files" "$OUTPUT_DIR/vcf_files" "$OUTPUT_DIR/read_mapping" "$OUTPUT_DIR/imputed_vcf_files" "$GMAP_DB_DIR"

if [[ "$DRY_RUN" -eq 0 ]]; then
    require_command phg
    require_command samtools
    require_command tabix
    if ! command -v "$GMAP_BUILD_CMD" >/dev/null 2>&1; then
        echo "Error: gmap_build not found at $GMAP_BUILD_CMD" >&2
        exit 1
    fi
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
    activate_conda
    phg --version
fi

echo "Using database root: $ROOT_DIR"
echo "Database name: $DATABASE_NAME"
echo "Reference fasta: $REFERENCE_FASTA"
echo "Reference gff: $REFERENCE_GFF"
echo "Assemblies: ${ASSEMBLY_SOURCES[*]}"

if [[ "$DRY_RUN" -eq 0 ]]; then
    run_cmd phg initdb --db-path "$VCF_DB_DIR/"
else
    echo "[DRY-RUN] phg initdb --db-path $VCF_DB_DIR/"
fi

REFERENCE_FASTA_FILE="$DATA_DIR/${REFERENCE_NAME}.fa"
REFERENCE_GFF_FILE="$DATA_DIR/${REFERENCE_NAME}.gff"

prepare_source_file "$REFERENCE_FASTA" "$REFERENCE_NAME" "$REFERENCE_FASTA_FILE"
prepare_source_file "$REFERENCE_GFF" "$REFERENCE_NAME" "$REFERENCE_GFF_FILE"

# Prepare the assemblies into the data directory
ASSEMBLY_FILES=()
for source in "${ASSEMBLY_SOURCES[@]}"; do
    sample_name="$(basename_without_ext "$source")"
    source_value="$source"
    if [[ "$source" == *"|"* ]]; then
        sample_name="${source%%|*}"
        source_value="${source#*|}"
    fi
    target_file="$DATA_DIR/${sample_name}.fa"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] would prepare assembly source: $source_value -> $target_file"
        ASSEMBLY_FILES+=("$target_file")
        continue
    fi
    prepare_source_file "$source_value" "$sample_name" "$target_file"
    ASSEMBLY_FILES+=("$target_file")
done

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would create reference ranges"
else
    run_cmd phg create-ranges --gff "$REFERENCE_GFF_FILE" --boundary gene -o "$OUTPUT_DIR/ref_ranges.bed" --reference-file "$REFERENCE_FASTA_FILE"
fi

ALIGNMENT_KEYFILE="$DATA_DIR/alignment_keyfile.txt"
rm -f "$ALIGNMENT_KEYFILE"
for assembly_file in "${ASSEMBLY_FILES[@]}"; do
    if [[ "$(basename "$assembly_file")" != "${REFERENCE_NAME}.fa" ]]; then
        echo "$assembly_file" >> "$ALIGNMENT_KEYFILE"
    fi
done

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would align assemblies using $ALIGNMENT_KEYFILE"
else
    run_cmd phg align-assemblies --gff "$REFERENCE_GFF_FILE" -o "$OUTPUT_DIR/alignment_files/" --total-threads "$THREADS" --in-parallel 2 --assembly-file-list "$ALIGNMENT_KEYFILE" --output-dir "$OUTPUT_DIR/alignment_files/" --reference-file "$REFERENCE_FASTA_FILE"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would build AGC archive"
else
    run_cmd phg agc-compress --db-path "$VCF_DB_DIR/" --fasta-list "$ALIGNMENT_KEYFILE" --reference-file "$REFERENCE_FASTA_FILE"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would create reference VCF"
else
    run_cmd phg create-ref-vcf --bed "$OUTPUT_DIR/ref_ranges.bed" --reference-file "$REFERENCE_FASTA_FILE" --reference-name "$REFERENCE_NAME" --db-path "$VCF_DB_DIR/"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would create MAF-based VCFs"
else
    for assembly_file in "${ASSEMBLY_FILES[@]}"; do
        assembly_name="$(basename_without_ext "$assembly_file")"
        if [[ "$assembly_name" == "$REFERENCE_NAME" ]]; then
            continue
        fi
        if [[ ! -f "$OUTPUT_DIR/vcf_files/${assembly_name}.h.vcf.gz" ]]; then
            run_cmd phg create-maf-vcf --bed "$OUTPUT_DIR/ref_ranges.bed" --reference-file "$REFERENCE_FASTA_FILE" -o "$OUTPUT_DIR/vcf_files/" --db-path "$VCF_DB_DIR/" --skip-metrics --maf-dir "$OUTPUT_DIR/alignment_files/"
            break
        fi
    done
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would index VCFs before loading"
    echo "[DRY-RUN] would load VCFs into the PHG database"
else
    index_vcf_files "$VCF_DB_DIR/hvcf_files" "$OUTPUT_DIR/vcf_files"
    run_cmd phg load-vcf --db-path "$VCF_DB_DIR/" --vcf-dir "$OUTPUT_DIR/vcf_files/" --threads "$THREADS"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would generate hapID files"
else
    if [[ ! -f "$OUTPUT_DIR/hapIDranges.tsv" ]]; then
        run_cmd phg sample-hapid-by-range --input-dir "$VCF_DB_DIR/hvcf_files/" --output-file "$OUTPUT_DIR/hapIDranges.tsv"
    fi
    if [[ ! -f "$OUTPUT_DIR/hapIDtable.tsv" ]]; then
        run_cmd phg hapid-sample-table --hvcf-dir "$VCF_DB_DIR/hvcf_files/" --output-file "$OUTPUT_DIR/hapIDtable.tsv"
    fi
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would build GMAP databases"
else
    for assembly_file in "${ASSEMBLY_FILES[@]}"; do
        assembly_name="$(basename_without_ext "$assembly_file")"
        genome_dir="$GMAP_DB_DIR/${assembly_name}"
        if [[ ! -f "$genome_dir/${assembly_name}.chromosome" ]]; then
            mkdir -p "$genome_dir"
            run_cmd "$GMAP_BUILD_CMD" -D "$GMAP_DB_DIR/" -d "$assembly_name" "$assembly_file"
        fi
    done
fi

HVCF2BED_URL="https://raw.githubusercontent.com/jsarriaa/PHGv2Tools/main/src/phgtools/modules/hvcf2bed.py"

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would download hvcf2bed.py from $HVCF2BED_URL to $ROOT_DIR/hvcf2bed.py"
    echo "[DRY-RUN] would create BED files from hVCFs"
else
    # download the official hvcf2bed module from GitHub if missing
    if [[ ! -f "$ROOT_DIR/hvcf2bed.py" ]]; then
        if command -v curl >/dev/null 2>&1; then
            run_cmd curl -fsSL "$HVCF2BED_URL" -o "$ROOT_DIR/hvcf2bed.py"
        elif command -v wget >/dev/null 2>&1; then
            run_cmd wget -qO "$ROOT_DIR/hvcf2bed.py" "$HVCF2BED_URL"
        else
            echo "Error: curl or wget required to download hvcf2bed.py from GitHub" >&2
            exit 1
        fi
    fi

    for assembly_file in "${ASSEMBLY_FILES[@]}"; do
        assembly_name="$(basename_without_ext "$assembly_file")"
        if [[ ! -f "$VCF_DB_DIR/hvcf_files/${assembly_name}.h.bed" ]]; then
            python3 "$ROOT_DIR/hvcf2bed.py" "$VCF_DB_DIR/hvcf_files/" "$assembly_name" --no-plot
        fi
    done
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] would index FASTA files"
else
    for assembly_file in "${ASSEMBLY_FILES[@]}"; do
        assembly_name="$(basename_without_ext "$assembly_file")"
        if [[ ! -f "$DATA_DIR/${assembly_name}.fa.fai" ]]; then
            samtools faidx "$assembly_file"
        fi
    done
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    if [[ -f "$ROOT_DIR/hvcf2bed.py" ]]; then
        echo "[DRY-RUN] hvcf2bed.py prepared at $ROOT_DIR/hvcf2bed.py"
    fi
else
    if [[ -f "$ROOT_DIR/hvcf2bed.py" ]]; then
        echo "Unified PHG database build completed for $DATABASE_NAME"
    fi
fi
