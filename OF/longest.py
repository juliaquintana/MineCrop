from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO

def parse_gff_attributes(attr_field: str) -> dict:
    """Parse a GFF3 attribute column into a dict, robustly (first '=' only, URL-decoded)."""
    attrs = {}
    for item in attr_field.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = unquote(value.strip())
        if key not in attrs:  # keep first occurrence
            attrs[key] = value
    return attrs

def strip_prefix(value: str, prefix: str) -> str:
    """Remove a known 'prefix:' like 'transcript:' or 'gene:' if present."""
    if value.startswith(prefix):
        return value[len(prefix):]
    return value

def normalize_id(x: str) -> str:
    """Strip whitespace, take first token, remove trailing version like '.1' if numeric."""
    x = x.strip().split()[0]
    if "." in x and x.rsplit(".", 1)[-1].isdigit():
        x = x.rsplit(".", 1)[0]
    return x

def main():
    script_dir = Path(__file__).resolve().parent
    gff_path = script_dir / "annotation.gff"
    faa_path = script_dir / "proteins.fa"
    out_path = script_dir / "longest_transcripts.pep"  # peptide FASTA output

    if not gff_path.exists():
        raise FileNotFoundError(f"GFF not found: {gff_path}")
    if not faa_path.exists():
        raise FileNotFoundError(f"FASTA not found: {faa_path}")

    # ---------- STEP 1: Build transcript -> gene map from mRNA rows ---------- #
    transcript_to_gene = {}
    total_gff = cds_gff = mrna_gff = 0

    with gff_path.open("r", encoding="utf-8", errors="replace") as gff:
        for line in gff:
            if not line or line.startswith("#"):
                continue
            total_gff += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype = fields[2]
            attrs = parse_gff_attributes(fields[8])

            if ftype == "mRNA":
                mrna_gff += 1
                tid = attrs.get("ID")
                parent_gene = attrs.get("Parent")
                if tid and parent_gene:
                    tid = strip_prefix(tid, "transcript:")
                    gid = strip_prefix(parent_gene, "gene:")
                    # Normalize (remove .1 suffixes)
                    tid = normalize_id(tid)
                    gid = normalize_id(gid)
                    transcript_to_gene[tid] = gid

    # ---------- STEP 2: Build protein -> gene map via CDS rows ---------- #
    protein_to_gene = {}
    with gff_path.open("r", encoding="utf-8", errors="replace") as gff:
        for line in gff:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype != "CDS":
                continue
            cds_gff += 1
            attrs = parse_gff_attributes(fields[8])

            pid = attrs.get("protein_id")
            parent_t = attrs.get("Parent")  # like 'transcript:ENSXJFT...'
            if not pid or not parent_t:
                continue

            pid = normalize_id(pid)
            tid = normalize_id(strip_prefix(parent_t, "transcript:"))

            gid = transcript_to_gene.get(tid)
            if gid:
                protein_to_gene[p_id := pid] = gid

    print(f"GFF rows: {total_gff}, mRNA rows: {mrna_gff}, CDS rows: {cds_gff}")
    print(f"Mappings: transcripts→genes={len(transcript_to_gene)}, proteins→genes={len(protein_to_gene)}")

    if not protein_to_gene:
        print("No protein→gene mappings found. Check that CDS rows have 'protein_id' and Parent=transcript:... "
              "and that mRNA rows map transcript IDs to gene IDs.")
        # Optionally: return early
        # return

    # ---------- STEP 3: Read proteins FASTA and keep longest per gene ---------- #
    longest = {}  # gene -> (record, length)
    fasta_records = matched = 0

    for record in SeqIO.parse(str(faa_path), "fasta"):
        fasta_records += 1
        # record.id is the first token of the header, e.g. 'ENSXJFP00005000015.1'
        raw_pid = record.id
        pid = normalize_id(raw_pid)

        # If FASTA 'gene:' and 'transcript:' are needed, they’re in record.description
        # but we rely on the protein_to_gene from the GFF for correctness.
        gid = protein_to_gene.get(pid)
        if not gid:
            # Some GFFs omit version on protein_id while FASTA has it; normalize_id above handles that.
            continue

        length = len(record.seq)
        matched += 1
        if gid not in longest or length > longest[gid][1]:
            longest[gid] = (record, length)

    print(f"FASTA records: {fasta_records}, matched to GFF proteins: {matched}, unique genes kept: {len(longest)}")

    # ---------- STEP 4: Write output ---------- #
    if longest:
        with out_path.open("w", encoding="utf-8") as out:
            for record, _ in longest.values():
                SeqIO.write(record, out, "fasta")
        print(f"Done! Wrote {len(longest)} longest sequences to {out_path}")
    else:
        print("No matches found—output would be empty. Verify ID normalization and mapping keys.")

if __name__ == "__main__":
    main()
