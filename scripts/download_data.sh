#!/bin/bash


# Download SRA data for Assignment 3
# Dataset: SRP126540
# 3 omnivore + 3 vegan samples


set -e  # stop on error

# WSL-native data directory (FAST)
DATA_DIR="$HOME/6110Assignment3/data/raw"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "Downloading SRA samples into $DATA_DIR"

# Number of threads (adjust if needed)
THREADS=6

# List of SRR IDs
SAMPLES=(
    SRR8146935
    SRR8146936
    SRR8146938
    SRR8146951
    SRR8146952
    SRR8146954
)

# Loop through samples
for srr in "${SAMPLES[@]}"
do
    echo "=============================="
    echo "Processing $srr..."
    echo "=============================="

    # Skip if already downloaded
    if [[ -f "${srr}_1.fastq.gz" ]]; then
        echo "$srr already downloaded, skipping..."
        continue
    fi

    # Download paired-end reads
    fasterq-dump "$srr" --split-files -e "$THREADS" -p

    # Compress FASTQ files
    gzip "${srr}_1.fastq"
    gzip "${srr}_2.fastq"

    echo "$srr completed."
done

echo "All downloads complete."