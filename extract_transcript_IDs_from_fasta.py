import re
import os

# === Configuration ===
# Use r"" for Windows paths to avoid escaping backslashes
fasta_path = r"C:\Users\julia.quintana\Baxter_Arabidopsis.fasta"

# Output name (will be saved in the same folder as the FASTA)
output_name = "transcript_ids.txt"

# === Regex: nuclear (1–5), mitochondria (M), chloroplast (C) ===
re_transcript = re.compile(r"AT(?:[1-5]|M|C)G\d{5}\.\d+")

def main():
    if not os.path.isfile(fasta_path):
        print(f"❌ FASTA not found: {fasta_path}")
        return

    transcript_ids = set()

    # Read FASTA and capture only header lines
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.lstrip().startswith(">"):
                transcript_ids.update(re_transcript.findall(line))

    sorted_ids = sorted(transcript_ids)

    # Build output path next to the FASTA
    fasta_dir = os.path.dirname(os.path.abspath(fasta_path)) or "."
    output_path = os.path.join(fasta_dir, output_name)

    # Write file
    with open(output_path, "w", encoding="utf-8") as out:
        for tid in sorted_ids:
            out.write(tid + "\n")

    print(f"✅ File written: {output_path}")
    print(f"✅ Total transcripts found: {len(sorted_ids)}")

    # Optional: show a few examples
    if sorted_ids:
        print("🔎 First 5 IDs:", ", ".join(sorted_ids[:5]))
    else:
        print("⚠️ No transcript IDs matched the regex. "
              "Check if headers contain ATxG#####.# or adjust the pattern.")

if __name__ == "__main__":
    main()