#!/usr/bin/env bash

# Exit immediately if an unexpected error occurs
set -e

# Expects path to graph folder ($1)
GRAPH_DIR="${1:-}"

# Check if a folder was provided and if it exists
if [ -z "$GRAPH_DIR" ] || [ ! -d "$GRAPH_DIR" ]; then
    echo "ERROR: No valid directory provided."
    echo "Usage: $0 <path_to_graph_folder>"
    exit 1
fi

# Remove trailing slash if present
GRAPH_DIR="${GRAPH_DIR%/}"

# Check if the assemblies.agc file exists inside the folder
AGC_FILE="${GRAPH_DIR}/assemblies.agc"
if [ ! -f "$AGC_FILE" ]; then
    echo "ERROR: Could not find 'assemblies.agc' at '${AGC_FILE}'."
    echo "Please ensure the path points to a valid graph folder containing a '.agc' file."
    echo "Usage: $0 <path_to_graph_folder>"
    exit 1
fi

# Ensure genomes.txt exists or can be generated if needed
GENOMES_FILE="${GRAPH_DIR}/genomes.txt"
if [ ! -f "$GENOMES_FILE" ]; then
    echo "Generating genomes list from assemblies.agc..."
    agc listset "$AGC_FILE" > "$GENOMES_FILE"
fi

mkdir -p "${GRAPH_DIR}/data"

while read g; do
    if [ ! -f "/gmap_db/${g}/${g}.chromosome" ]; then
        echo "formatting ${g}"
        agc getset -o "${GRAPH_DIR}/data/${g}.fa" "$AGC_FILE" "$g"
        /barleygraph/gmap-2013-08-31/local/bin/gmap_build -D /gmap_db/ -d "${g}" "${GRAPH_DIR}/data/${g}.fa" &> "/gmap_db/${g}.log"
        
        # Remove ONLY the specific fasta generated for this loop iteration
        rm -f "${GRAPH_DIR}/data/${g}.fa"
    fi
done < "$GENOMES_FILE"