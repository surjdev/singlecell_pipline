#!/usr/bin/env python3
"""
Generate a Nextflow-compatible samplesheet from docs/filenames.txt
Usage:
    python scripts/generate_samplesheet.py [filenames_txt] [base_dir] [out_csv]
"""
import csv
import sys
from pathlib import Path

def main():
    filenames_file = Path(sys.argv[1] if len(sys.argv) > 1 else "docs/filenames.txt")
    base_dir = sys.argv[2] if len(sys.argv) > 2 else "/home/koraop/nob_dir/sc_transcriptomics/fastq_out"
    out_csv = Path(sys.argv[3] if len(sys.argv) > 3 else "samplesheet.csv")

    if not filenames_file.exists():
        print(f"Error: {filenames_file} not found!", file=sys.stderr)
        sys.exit(1)

    with open(filenames_file, "r") as f:
        tokens = f.read().split()

    # Filter only fastq.gz files
    fastq_files = [t for t in tokens if t.endswith(".fastq.gz") or t.endswith(".fq.gz")]
    fastq_files.sort()

    rows = []
    for fname in fastq_files:
        # Extract sample ID: e.g. HP1525301T2D_L17.fastq.gz -> HP1525301T2D_L17
        sample_id = fname.replace(".fastq.gz", "").replace(".fq.gz", "")
        full_path = f"{base_dir.rstrip('/')}/{fname}"
        rows.append({
            "sample": sample_id,
            "fastq_1": full_path,
            "fastq_2": ""  # Single-End: fastq_2 is empty
        })

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "fastq_1", "fastq_2"])
        writer.writeheader()
        writer.writerows(rows)

    print(f" Generated {len(rows)} single-cell samples into: {out_csv}")
    print(f"   Base path: {base_dir}")
    print("   Format: Single-End (fastq_2 is blank)")

if __name__ == "__main__":
    main()
