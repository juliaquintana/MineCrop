from Bio import SeqIO
import os

def extract_sequences(fasta_file, id_list_file, output_file):

    # Check ID file exists
    if not os.path.isfile(id_list_file):
        print(f"❌ ERROR: ID list file not found:\n{id_list_file}")
        return

    # Check FASTA exists
    if not os.path.isfile(fasta_file):
        print(f"❌ ERROR: FASTA file not found:\n{fasta_file}")
        return

    # Load partial IDs
    with open(id_list_file, 'r') as f:
        partial_ids = set(line.strip() for line in f)

    print(f"Loaded {len(partial_ids)} IDs")

    found = 0
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if any(pid in record.id for pid in partial_ids):
                SeqIO.write(record, out_f, "fasta")
                found += 1

    print(f"✅ Extraction complete. Sequences found: {found}")
    print(f"📄 Output written to: {output_file}")


# Example usage
if __name__ == "__main__":
    extract_sequences(
        r"C:\Users\julia.quintana\input.fasta",
        r"C:\Users\julia.quintana\ids.txt",
        r"C:\Users\julia.quintana\output.fasta"
    )


