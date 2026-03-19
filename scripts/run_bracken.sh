#!/bin/bash

# ===============================
# Bracken abundance estimation
# ===============================

set -e

# Paths
DB_DIR="$HOME/kraken_db"
KRAKEN_DIR="$HOME/kraken_results"
OUT_DIR="$HOME/bracken_results"

mkdir -p "$OUT_DIR"

# Taxonomic level
LEVEL="S"   # S = species (you can also do G = genus)

# Read length
READ_LEN=150

# Samples
SAMPLES=(
    SRR8146935
    SRR8146936
    SRR8146938
    SRR8146951
    SRR8146952
    SRR8146954
)

echo "Starting Bracken..."

for srr in "${SAMPLES[@]}"
do
    echo "=============================="
    echo "Processing $srr..."
    echo "=============================="

    REPORT="${KRAKEN_DIR}/${srr}_report.txt"
    OUTPUT="${OUT_DIR}/${srr}_bracken.txt"

    if [[ -f "$OUTPUT" ]]; then
        echo "$srr already processed, skipping..."
        continue
    fi

    bracken \
        -d "$DB_DIR" \
        -i "$REPORT" \
        -o "$OUTPUT" \
        -r "$READ_LEN" \
        -l "$LEVEL"

    echo "$srr completed."
done

echo "All Bracken runs complete."