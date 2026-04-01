#!/usr/bin/env python3

# Input files (adjust paths if needed)
ids_file = "/Users/julia/Baxter_Araport11_pep_IDs.txt"   # one ID per line
og_file = "/Users/julia/orthogroups_with_arabidopsis"                 # OG table
out_file = "/Users/julia/filtered_OG.tsv"               # output file

# Load IDs (one per line)
with open(ids_file) as f:
    ids = {line.strip() for line in f if line.strip()}

print(f"Loaded {len(ids)} IDs from {ids_file}")

# Scan orthogroup file, keep lines containing any of the IDs
with open(og_file) as fin, open(out_file, "w") as fout:
    for line in fin:
        if any(id_ in line for id_ in ids):
            fout.write(line)

print(f"Done. Matching lines written to {out_file}")
