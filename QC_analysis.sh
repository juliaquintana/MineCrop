#!/bin/bash

# Description: Download FASTQ files from SRA and run QC with fastp
# Author: Julia Quintana
# Date: 2025-06-24

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRA_LIST="$SCRIPT_DIR/SRR_Acc_List.txt"
RAW_DIR="$SCRIPT_DIR/raw_fastq"
QC_DIR="$SCRIPT_DIR/qc_reports"

# Create output directories
mkdir -p "$RAW_DIR" "$QC_DIR"

# Loop through each accession number
while IFS= read -r accession; do
    echo "🔽 Downloading $accession..."
    prefetch "$accession" && \
    fasterq-dump "$accession" -O "$RAW_DIR" -t "$RAW_DIR"

    echo "🔬 Running fastp QC for $accession..."
    fastp \
        -i "$RAW_DIR/${accession}_1.fastq" \
        -I "$RAW_DIR/${accession}_2.fastq" \
        -o "$QC_DIR/${accession}_1.clean.fastq" \
        -O "$QC_DIR/${accession}_2.clean.fastq" \
        -h "$QC_DIR/${accession}_fastp.html" \
        -j "$QC_DIR/${accession}_fastp.json" \
        --thread 4

# Upload to Google Cloud Storage
        gsutil -m cp -r "$QC_DIR/${accession}"* gs://$QC_DIR

# Clear the cache
        cache-mgr --clear
        echo "Cleared cache for $accession."
    else
        echo "❌ Failed to prefetch $accession"
    fi
    echo "✅ Done with $accession"
done < "$SRA_LIST"
