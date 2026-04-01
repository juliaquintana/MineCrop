
<p align="center">
  <img src="https://github.com/user-attachments/assets/159906b0-8cd0-4c9a-999f-2a18c5fbc977" width="500" />
</p>

# MineCrop

**MineCrop** is a research project funded by the **Comunidad de Madrid** and hosted by **Universidad Rey Juan Carlos (URJC)**.  
The project is firmly committed to **open science**, **open data**, and the public sharing of computational resources developed throughout its execution.

## Project Overview

The MineCrop project aims to enhance our understanding of **crop domestication** by identifying **micronutrient metabolism pathways** that may have been **lost or altered during the domestication of traditional crop varieties**.

By combining comparative genomics, transcriptomics, and functional annotation, the project seeks to identify major genetic changes associated with nutritional traits and domestication processes.

All scripts developed during the project are made **freely available** in this repository and are designed to automate specific data handling, processing, and analysis workflows commonly used in genomics and evolutionary biology.

## License

This repository is licensed under the MIT License.

---

## Repository Contents

This repository currently includes:

1. **Shell scripts** for downloading, processing, and quality‑controlling SRA sequencing data  
2. **Python scripts** for orthogroup analysis and post‑processing after running **OrthoFinder**

---

## 1. Shell Scripts for Downloading and Processing SRA Data

### `sra_fastq_qc_pipeline.sh`

#### Purpose
Downloads sequencing data from the **NCBI Sequence Read Archive (SRA)**, converts it to FASTQ format, performs quality control and trimming using **fastp**, and uploads cleaned reads and QC reports to **Google Cloud Storage**.

#### Description
This script processes SRA accessions **sequentially** from a provided list. For each accession, it:

- Downloads raw SRA files using `prefetch`
- Converts SRA files to paired‑end FASTQ using `fasterq-dump`
- Runs **fastp** for quality control, adapter detection, and read trimming
- Generates compressed cleaned FASTQ files and fastp HTML/JSON reports
- Uploads QC outputs to a specified Google Cloud Storage bucket
- Cleans up intermediate FASTQ files and clears the SRA cache

Temporary files are stored in `/tmp` to improve I/O performance and reduce filesystem clutter.

#### Requirements
- Bash (Linux or macOS)
- SRA Toolkit (`prefetch`, `fasterq-dump`, `cache-mgr`)
- fastp
- Google Cloud SDK (`gsutil`)
- Internet access
- Sufficient disk space in `/tmp`

#### Input
- A plain text file containing one SRA accession per line  
  (`SRR_Acc_List.txt`)

#### Output
- Cleaned paired‑end FASTQ files (`*.clean.fastq.gz`)
- fastp QC reports (`*.html`, `*.json`)
- Uploaded results in a Google Cloud Storage bucket

#### Typical Use Cases
- Quality control of small to medium numbers of SRA datasets
- Simple and transparent SRA‑to‑FASTQ workflows
- Environments where parallel execution is not required
- Step‑by‑step debugging or manual inspection of individual samples

---

### `sra_fastq_qc_pipeline_parallel.sh`

#### Purpose
Performs **parallelized** downloading, quality control, and trimming of SRA datasets with enhanced robustness, logging, and error handling.

#### Description
This script is a **production‑ready**, scalable version of the SRA QC pipeline that processes multiple accessions in parallel using **GNU parallel**. It includes:

- Safe execution (`set -euo pipefail`)
- Timestamped logging for traceability
- Automatic retry logic for network‑dependent commands
- Success and failure tracking
- Centralized cleanup after all jobs complete

For each SRA accession, the pipeline:
- Downloads the SRA record using `prefetch` with retries
- Converts it to paired‑end FASTQ format
- Performs QC and trimming with **fastp**
- Uploads cleaned reads and QC reports to Google Cloud Storage
- Records the accession as successful or failed

#### Requirements
- Bash (Linux recommended)
- GNU parallel
- SRA Toolkit (`prefetch`, `fasterq-dump`, `cache-mgr`)
- fastp
- Google Cloud SDK (`gsutil`)
- Internet access
- Multiple CPU cores
- Sufficient temporary storage in `/tmp`

#### Input
- A text file with one SRA accession per line  
  (`SRR_Acc_List.txt`)

#### Output
- Cleaned paired‑end FASTQ files (`*.clean.fastq.gz`)
- fastp QC reports (`*.html`, `*.json`)
- Uploaded results in Google Cloud Storage
- Summary files:
  - `successful_accessions.txt`
  - `failed_accessions.txt`

#### Typical Use Cases
- High‑throughput SRA data processing
- Large RNA‑seq or genomic projects
- Cloud or HPC environments
- Automated pipelines requiring robustness, scalability, and reporting

---

## 2. Python Scripts for Orthogroup and Comparative Genomics Analysis

These scripts are primarily used **after running OrthoFinder**, to subset, filter, and summarize gene and protein data.

---

### `extract_fasta.py`

#### Purpose
Extracts sequences from a FASTA file whose identifiers match a list of provided IDs. Partial matches are supported.

#### Description
This script reads:
- A FASTA file
- A list of IDs (one per line)

Any sequence whose record ID contains one of the listed identifiers is extracted and written to a new FASTA file. This is useful for subsetting large reference FASTA files based on genes, transcripts, or proteins of interest.

#### Requirements
- Python 3
- Biopython

#### Typical Use Cases
- Extracting specific transcripts or proteins from reference FASTA files
- Creating reduced FASTA datasets for downstream analyses

---

### `extract_transcript_IDs_from_fasta.py`

#### Purpose
Extracts transcript IDs from FASTA headers using a regular expression.

#### Description
The script scans FASTA header lines and extracts **Arabidopsis‑style transcript IDs**
(e.g. `AT1G01010.1`, including nuclear, mitochondrial, and chloroplast IDs).  
Extracted IDs are written to a text file, one per line.

#### Typical Use Cases
- Generating transcript ID lists from reference FASTA files
- Preparing ID files for sequence extraction or orthogroup filtering

#### Notes
- The transcript ID pattern is defined using a regular expression and can be modified if needed.

---

### `filter_orthogroups.py`

#### Purpose
Filters an orthogroup table to retain only rows containing specified sequence IDs.

#### Description
This script reads:
- A list of sequence IDs
- An orthogroup table (e.g. OrthoFinder output)

Any line containing one or more of the provided IDs is retained and written to a new output file.

#### Typical Use Cases
- Identifying orthogroups containing genes of interest
- Reducing large orthogroup tables to targeted subsets

#### Notes
- Input and output file paths are currently hardcoded and should be adjusted before execution.

---

### `longest.py`

#### Purpose
Selects the **longest protein isoform per gene** using GFF annotations.

#### Description
This script combines information from:
- A GFF annotation file
- A protein FASTA file

For genes with multiple protein isoforms, the script retains only the **longest protein sequence**.

#### Workflow Summary
- Build transcript → gene mappings from mRNA features in the GFF
- Build protein → gene mappings from CDS features
- Parse protein FASTA sequences
- Retain the longest protein per gene

#### Requirements
- Python 3
- Biopython

#### Typical Use Cases
- Generating representative, non‑redundant protein datasets
- Preparing input FASTA files for comparative genomics analyses

---

### `og_consensus.py`

#### Purpose
Generates **consensus functional annotations** for orthogroups based on EggNOG‑mapper results.

#### Description
The script reads multiple `*.emapper.annotations` files (one per orthogroup) and summarizes functional annotations across all member genes. Consensus annotations are derived using frequency‑based thresholds.

#### Outputs
- A TSV file with one row per orthogroup containing:
  - Functional description
  - Preferred gene name
  - GO terms
  - KEGG pathways
  - EC numbers
  - PFAM domains
  - COG categories
- Optional TSV file with per‑gene annotation details

#### Requirements
- Python 3
- pandas

#### Typical Use Cases
- Functional characterization of orthogroups
- Comparative genomics and enrichment analyses

---

## General Notes

- All Python scripts are written in **Python 3**
- Input paths, parameters, and thresholds should be reviewed and adjusted before execution
- Scripts are intended for **Linux/macOS** or other Unix‑like environments
- Contributions, reuse, and adaptation are encouraged in the spirit of **open science**

---

## Citation & Reuse

If you use this repository or its scripts in published work, please cite the MineCrop project and acknowledge **Comunidad de Madrid** and **Universidad Rey Juan Carlos**.

---

``
