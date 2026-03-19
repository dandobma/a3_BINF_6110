#!/bin/bash


# Kraken2 classification script


set -e

# Paths
DATA_DIR="$HOME/6110Assignment3/data"
DB_DIR="$HOME/kraken_db"
OUT_DIR="$HOME/kraken_results"

mkdir -p "$OUT_DIR"

# Threads
THREADS=6

# Samples
SAMPLES=(
    SRR8146935
    SRR8146936
    SRR8146938
    SRR8146951
    SRR8146952
    SRR8146954
)

echo "Starting Kraken2 classification..."

for srr in "${SAMPLES[@]}"
do
    echo "=============================="
    echo "Processing $srr..."
    echo "=============================="

    # Input files
    R1="${DATA_DIR}/${srr}_1.fastq.gz"
    R2="${DATA_DIR}/${srr}_2.fastq.gz"

    # Output files
    REPORT="${OUT_DIR}/${srr}_report.txt"
    OUTPUT="${OUT_DIR}/${srr}_output.txt"

    # Skip if already done
    if [[ -f "$REPORT" ]]; then
        echo "$srr already processed, skipping..."
        continue
    fi

    kraken2 \
        --db "$DB_DIR" \
        --threads "$THREADS" \
        --paired \
        --gzip-compressed \
        --report "$REPORT" \
        --output "$OUTPUT" \
        "$R1" "$R2"

    echo "$srr completed."
done

echo "All Kraken2 runs complete."