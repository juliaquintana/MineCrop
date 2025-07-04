#!/bin/bash
# Description: Download FASTQ files from SRA and run QC with fastp, including trimming
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
    
    if prefetch "$accession"; then
        fasterq-dump "$accession" -O "$RAW_DIR" -t "$RAW_DIR"
        echo "🔬 Running fastp QC for $accession..."
        # Run fastp QC analysis for PE reads. For SE reads, use only -i and -o and remove _1 suffix. --trimm is usually done after a first QC check
        fastp \
            -i "$RAW_DIR/${accession}_1.fastq" \
            -I "$RAW_DIR/${accession}_2.fastq" \
            -o "$QC_DIR/${accession}_1.clean.fastq" \
            -O "$QC_DIR/${accession}_2.clean.fastq" \
            -h "$QC_DIR/${accession}_fastp.html" \
            -j "$QC_DIR/${accession}_fastp.json" \
            --detect_adapter_for_pe \
            --trim_front1 15 \
            --trim_front2 15 \
            --thread 4

        # Compress cleaned FASTQ files
        echo "📦 Compressing cleaned FASTQ files..."
        gzip "$QC_DIR/${accession}_1.clean.fastq"
        gzip "$QC_DIR/${accession}_2.clean.fastq"
        
        # Upload to Google Cloud Storage
        echo "☁️ Uploading QC results to Google Cloud..."
        gsutil -m cp -r "$QC_DIR/${accession}_fastp.html" "$QC_DIR/${accession}_1.clean.fastq.gz" "$QC_DIR/${accession}_2.clean.fastq.gz" gs://$QC_DIR_in_Google_Cloud/QC

        # Clear the cache
        echo "🧹 Clearing cache..."
        cache-mgr --clear
        echo "Cleared cache for $accession."
    else
        echo "❌ Failed to prefetch $accession"
    fi

    echo "✅ Done with $accession"
done < "$SRA_LIST"
