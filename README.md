**MineCrop**

<img src="https://github.com/user-attachments/assets/159906b0-8cd0-4c9a-999f-2a18c5fbc977" width="500">

This is the repository for the project MineCrop, funded by Comunidad de Madrid and hosted by Universidad Rey Juan Carlos. This project has a commitement with open data and shared resources.

This project aims to enhance our understanding of crop domestication by identifying micronutrient metabolism pathways that may have been lost during domestication of traditional varieties.

All the scripts developed during the project will be uploaded here and made publicly available. Each script is designed to automate specific tasks related to data handling and processing workflows.

2. Python scripts used for sequence processing and orthogroup analysis
**extract_fasta.py**
Purpose:
Extracts sequences from a FASTA file whose identifiers match a list of provided IDs (partial matches are supported).
Description:
This script reads a FASTA file and an ID list (one ID per line). Any sequence whose record ID contains one of the listed IDs is extracted and written to a new FASTA file. It is useful for subsetting large FASTA files based on transcript, gene, or protein identifiers.

Requirements:
Python 3
Biopython

Typical use cases:
Extracting specific transcripts or proteins from a reference FASTA
Creating reduced FASTA files for downstream analysis

**extract_transcript_IDs_from_fasta.py**
Purpose:
Extracts transcript IDs from FASTA headers using a regular expression.
Description:
This script scans FASTA header lines and extracts Arabidopsis-style transcript IDs (e.g. AT1G01010.1, including nuclear, mitochondrial, and chloroplast IDs). The extracted IDs are written to a text file, one per line.

Typical use cases:
Generating transcript ID lists from reference FASTA files
Preparing input ID files for sequence extraction or orthogroup filtering

Notes:
The transcript ID pattern is defined using a regular expression and can be modified if needed.

**filter_orthogroups.py**
Purpose:
Filters an orthogroup table to retain only rows containing specified sequence IDs.
Description:
The script reads a list of IDs and scans an orthogroup table file. Any line containing one or more of the provided IDs is written to a new output file. This allows selection of orthogroups associated with specific genes or proteins.

Typical use cases:
Identifying orthogroups containing genes of interest
Reducing large orthogroup tables to a targeted subset

Notes:
Input and output file paths are currently hardcoded and should be adjusted before running.

**longest.py**
Purpose:
Selects the longest protein sequence per gene based on GFF annotation.

Description:
This script parses a GFF annotation file and a protein FASTA file to determine which protein corresponds to each gene. For genes with multiple protein isoforms, only the longest protein sequence is retained and written to a new FASTA file.

Workflow summary:
Build transcript → gene mappings from mRNA features in the GFF
Build protein → gene mappings from CDS features
Parse protein FASTA
Keep the longest protein sequence per gene

Requirements:
Python 3
Biopython

Typical use cases:
Generating representative protein sets
Preparing non-redundant protein FASTA files for comparative analyses

**og_consensus.py**
Purpose:
Generates consensus functional annotations for orthogroups based on EggNOG-mapper results.

Description:
This script reads multiple *.emapper.annotations files (one per orthogroup), summarizes functional annotations across all member genes, and produces a consensus annotation per orthogroup. It reports consensus descriptions, preferred names, GO terms, KEGG annotations, EC numbers, PFAM domains, and COG categories based on frequency thresholds.

Outputs:
A TSV file with one row per orthogroup containing consensus annotations
An optional TSV file containing per-gene annotation details for all orthogroups

Requirements:
Python 3
pandas

Typical use cases:
Functional characterization of orthogroups
Downstream comparative genomics or enrichment analyses

*General Notes
All scripts are written in Python 3.
Check and adjust input paths, parameters, and thresholds before execution.
Scripts are intended for use in Linux/macOS or other Unix-like environments.*
