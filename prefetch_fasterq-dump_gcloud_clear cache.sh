#!/bin/bash

# Description: This script automatically downloads NCBI SRA reads and stores them in Google Cloud to avoid disk space issues. This script is meant to be used after installing SRAtoolkit in your computer.
# Author: Julia Quintana
# Date: 2025-02-03
# Version: 1.0

# File containing the list of SRA accession numbers in the same folder as the script
SRA_LIST="/path/to/SRR_Acc_List.txt"

# Loop through each accession number in the list. For fasterq-dump, specify a folder with enough space for temporary .temp and .fastq files individual downloads. After uploading to Google Storage, the files will be removed.
while IFS= read -r accession; do
    echo "Processing $accession"
    if prefetch.3.2.0 "$accession"; then
        fasterq-dump.3.2.0 "$accession" --outdir "/path/to/sra" -t "/path/to/sra"
        # Upload to Google Cloud Storage
        gsutil -m cp -r "/path/to/sra/${accession}*" gs://srastorage/
        # Clear the cache after uploading to Google Storage
        cache-mgr --clear
        echo "Cleared cache for $accession."
    else
        echo "Failed to prefetch $accession"
    fi
done < "$SRA_LIST"

# The end
