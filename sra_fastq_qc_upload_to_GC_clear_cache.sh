#!/bin/bash
# Description: Download FASTQ files from SRA and run QC with fastp, including trimming
# Author: Julia Quintana
# Date: 2025-06-24
# Temporary storage: Uses /tmp for intermediate files to reduce clutter and improve I/O speed.
# Compression: Ensures output is compressed with .gz to save space.
# Whenever possible, use more threads to increase spped

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRA_LIST="$SCRIPT_DIR/SRR_Acc_List.txt"
# Temporary storage for raw files: Uses /tmp for intermediate files to reduce clutter and improve I/O speed. After QC, those are no longer needed.
RAW_DIR="/tmp/sra_raw"
QC_DIR="$SCRIPT_DIR/qc_reports"

# Create output directories
mkdir -p "$RAW_DIR" "$QC_DIR"

# Loop through each accession number
while IFS= read -r accession; do
    echo "🔽 Downloading $accession..."
    
    if prefetch "$accession"; then
        fasterq-dump "$accession" -O "$RAW_DIR" -t "$RAW_DIR"
        echo "🔬 Running fastp QC for $accession..."
        # Run fastp QC
        fastp \
            -i "$RAW_DIR/${accession}_1.fastq" \
            -I "$RAW_DIR/${accession}_2.fastq" \
            -o "$QC_DIR/${accession}_1.clean.fastq.gz" \
            -O "$QC_DIR/${accession}_2.clean.fastq.gz" \
            -h "$QC_DIR/${accession}_fastp.html" \
            -j "$QC_DIR/${accession}_fastp.json" \
            --detect_adapter_for_pe \
            --trim_front1 15 \
            --trim_front2 15 \
            --thread 4
        
        # Upload to Google Cloud Storage. You should use the -r (recursive) flag with gsutil cp only when you're copying directories, not individual files.
        echo "☁️ Uploading QC results to Google Cloud..."
        gsutil -m cp "$QC_DIR/${accession}_fastp.html" \
             "$QC_DIR/${accession}_1.clean.fastq.gz" \
             "$QC_DIR/${accession}_2.clean.fastq.gz" \
             gs://srastorage/Triticum_shoot

        # Clear the cache
                echo "🧹 Cleaning up..."
        cache-mgr --clear
        rm -f "$RAW_DIR/${accession}_1.fastq" "$RAW_DIR/${accession}_2.fastq"
    else
        echo "❌ Failed to prefetch $accession"
    fi

    echo "✅ Done with $accession"
done < "$SRA_LIST"
