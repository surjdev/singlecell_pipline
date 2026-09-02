"""
Single-Cell FASTQ & Whitelist Synthetic Data Generator.
Generates realistic 10x Chromium 3' v3 scRNA-seq FASTQ paired files (R1 & R2)
incorporating cell barcodes, UMIs, cDNA sequences, adapter artifacts,
poly-A/poly-G tails, low-quality ends, and ambient noise.
"""

import gzip
import random
import argparse
from pathlib import Path
from typing import List, Tuple
import numpy as np
from rich.progress import track
from rich.console import Console

console = Console()

BASES = ['A', 'C', 'G', 'T']
TSO_ADAPTER = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
POLY_A = "A" * 30
POLY_G = "G" * 30

MOCK_GENES = [
    ("GAPDH", "GTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTC"),
    ("ACTB", "GCTGTGCTATGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGG"),
    ("MALAT1", "AGGCGTTTGGGGAAAGCGTGGTGTTAGTGTTAGTGTTCCTCGGTCCACGTGCGTGTCTTTGTTTTAGTCTTTCCTTGTCCTTTTGTCCTCT"),
    ("CD3D", "ATGGAACATAGCACGTTTCTCTCTGGCCTGGTACTGGCTACCCTTCTCTCGCAAGTGAGCCCCTTCAAGATACCTATAGAGGAACTTGAG"),
    ("CD19", "ATGCTCCCTGGGGCCATGGGCAGCAGTGCCACCGTAGGGGTCTTCCTGCTACTGCTACTGCTACTGCCCGGCCCGTCACCGGCCGCCGCC"),
    ("MS4A1", "ATGACAACACCCAGAAATTCAGTAAATGGGACTTTCCCGGCAGAGCCAATGAAAGGCCCTATTGCTATGCAATCTGGTCCAAAACCACTC"),
    ("IL7R", "ATGATGATGATCATGTTGCTCAGTGTTCTTCTTGGTGGGGGTATGTACTTGCTGCCTCAGACTCCCGTTCTCCGGGAAGAGATGACCATG"),
    ("NKG7", "ATGGAGAGCGTCATGGTGCTGCTGCTGGTGGCCATCATCGTGCTGGCCGTGGCGCTCCCGGTGCCGACGCAGGGCTCGGTGCAGGAGGAG"),
    ("CST3", "ATGGCCCGGCCCCTGCGTGCCCTGCTCGCCCTCCTGGCCCTGGCCGCCCTGGCCATCAGCGCCAGCGAGGACAAGTTTCTCCGGCGTCGG"),
    ("PPBP", "ATGAGCTCCACCGCAGTCACCTTCTTCCTGCTCCTGCTGCTGACGCTGCAGCTGGGGCCCGCCGCCCCTGCCCCCCGTTTCCCACTTCGC")
]

def generate_random_dna(length: int) -> str:
    """Generate a random DNA string of given length."""
    return "".join(random.choices(BASES, k=length))

def generate_quality_string(length: int, mean_phred: int = 35, std_phred: int = 3, degrade_tail: bool = False) -> str:
    """Generate Phred+33 ASCII quality string with optional 3' degradation."""
    quals = np.random.normal(mean_phred, std_phred, length)
    if degrade_tail and length > 20:
        decay = np.linspace(0, 15, length - 20)
        quals[20:] -= decay
    quals = np.clip(quals, 10, 40).astype(int)
    return "".join(chr(q + 33) for q in quals)

def mutate_barcode(cb: str, max_dist: int = 1) -> str:
    """Introduce 1-base substitution error into a barcode."""
    cb_list = list(cb)
    pos = random.randint(0, len(cb) - 1)
    orig_base = cb_list[pos]
    cb_list[pos] = random.choice([b for b in BASES if b != orig_base])
    return "".join(cb_list)

def generate_whitelist(num_barcodes: int = 1000, cb_len: int = 16) -> List[str]:
    """Generate a realistic set of unique cell barcodes."""
    barcodes = set()
    while len(barcodes) < num_barcodes:
        barcodes.add(generate_random_dna(cb_len))
    return sorted(list(barcodes))

def generate_single_cell_dataset(
    outdir: str = "data/raw",
    whitelist_file: str = "data/whitelist/737K-august-2016.txt",
    sample_id: str = "scRNA_sample_01",
    num_cells: int = 100,
    reads_per_cell: int = 200,
    num_ambient_cells: int = 300,
    cb_len: int = 16,
    umi_len: int = 12,
    cdna_len: int = 91,
    seed: int = 42
) -> Tuple[Path, Path, Path]:
    """Generate paired 10x FASTQ.gz files and whitelist."""
    random.seed(seed)
    np.random.seed(seed)
    
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)
    wl_path = Path(whitelist_file)
    wl_path.parent.mkdir(parents=True, exist_ok=True)
    
    r1_file = out_path / f"{sample_id}_R1.fastq.gz"
    r2_file = out_path / f"{sample_id}_R2.fastq.gz"
    
    console.print(f"[bold cyan]Generating mock 10x scRNA-seq dataset...[/bold cyan]")
    
    # 1. Generate full whitelist (including uncaptured background)
    total_wl_size = max(1000, num_cells * 10)
    whitelist = generate_whitelist(total_wl_size, cb_len)
    with open(wl_path, "w") as f:
        for bc in whitelist:
            f.write(f"{bc}\n")
    console.print(f"  ✔ Generated whitelist: [green]{wl_path}[/green] ({len(whitelist):,} barcodes)")
    
    # Select true cell barcodes vs ambient background
    true_cells = whitelist[:num_cells]
    ambient_barcodes = whitelist[num_cells:num_cells + num_ambient_cells]
    
    # Cell read depths follow log-normal distribution
    cell_depths = np.random.lognormal(mean=np.log(reads_per_cell), sigma=0.4, size=num_cells).astype(int)
    ambient_depths = np.random.geometric(p=0.2, size=num_ambient_cells) # mostly 1-5 reads
    
    total_expected_reads = sum(cell_depths) + sum(ambient_depths)
    console.print(f"  ✔ Simulating {num_cells} real cells + {num_ambient_cells} ambient droplets ({total_expected_reads:,} total reads)")
    
    read_idx = 0
    with gzip.open(r1_file, "wt") as f1, gzip.open(r2_file, "wt") as f2:
        # Simulate real cells
        for cell_i, (cb, depth) in enumerate(zip(true_cells, cell_depths)):
            for _ in range(depth):
                read_idx += 1
                read_name = f"@{sample_id}:{read_idx}:FLOWCELL:1:1101:{random.randint(1000,9999)}:{random.randint(1000,9999)}"
                
                # Barcode error simulation (5% 1-bp mismatch, 1% invalid random)
                p = random.random()
                if p < 0.05:
                    r1_cb = mutate_barcode(cb)
                elif p < 0.06:
                    r1_cb = generate_random_dna(cb_len)
                else:
                    r1_cb = cb
                    
                # UMI (12bp)
                r1_umi = generate_random_dna(umi_len)
                r1_seq = r1_cb + r1_umi
                r1_qual = generate_quality_string(cb_len + umi_len, mean_phred=37, std_phred=2)
                
                # cDNA simulation (R2)
                gene_name, gene_seq = random.choice(MOCK_GENES)
                # Introduce biological transcript variability & artifacts
                art_prob = random.random()
                if art_prob < 0.08:
                    # Adapter / TSO read-through contamination
                    r2_seq = gene_seq[:40] + TSO_ADAPTER + generate_random_dna(21)
                elif art_prob < 0.15:
                    # Poly-A tail at 3' end
                    r2_seq = gene_seq[:60] + POLY_A
                elif art_prob < 0.20:
                    # Poly-G NextSeq 2-color dark cycle artifact
                    r2_seq = gene_seq[:50] + POLY_G
                else:
                    # Normal transcript sequence
                    offset = random.randint(0, max(0, len(gene_seq) - cdna_len))
                    sub_seq = gene_seq[offset:offset + cdna_len]
                    if len(sub_seq) < cdna_len:
                        sub_seq += generate_random_dna(cdna_len - len(sub_seq))
                    r2_seq = sub_seq
                    
                r2_seq = r2_seq[:cdna_len]
                # 3' quality degradation on R2
                r2_qual = generate_quality_string(len(r2_seq), mean_phred=34, std_phred=4, degrade_tail=True)
                
                # Write R1
                f1.write(f"{read_name} 1:N:0:0\n{r1_seq}\n+\n{r1_qual}\n")
                # Write R2
                f2.write(f"{read_name} 2:N:0:0\n{r2_seq}\n+\n{r2_qual}\n")
                
        # Simulate ambient droplets
        for cb, depth in zip(ambient_barcodes, ambient_depths):
            for _ in range(depth):
                read_idx += 1
                read_name = f"@{sample_id}:{read_idx}:FLOWCELL:1:1101:{random.randint(1000,9999)}:{random.randint(1000,9999)}"
                r1_cb = cb
                r1_umi = generate_random_dna(umi_len)
                r1_seq = r1_cb + r1_umi
                r1_qual = generate_quality_string(cb_len + umi_len, mean_phred=32, std_phred=4)
                
                gene_name, gene_seq = random.choice(MOCK_GENES[:3]) # Ambient RNA dominated by high abundance genes
                r2_seq = gene_seq[:cdna_len]
                if len(r2_seq) < cdna_len:
                    r2_seq += generate_random_dna(cdna_len - len(r2_seq))
                r2_qual = generate_quality_string(len(r2_seq), mean_phred=30, std_phred=5, degrade_tail=True)
                
                f1.write(f"{read_name} 1:N:0:0\n{r1_seq}\n+\n{r1_qual}\n")
                f2.write(f"{read_name} 2:N:0:0\n{r2_seq}\n+\n{r2_qual}\n")

    console.print(f"[bold green]✔ Done![/bold green] Generated {read_idx:,} paired reads in:")
    console.print(f"  • R1: [cyan]{r1_file}[/cyan]")
    console.print(f"  • R2: [cyan]{r2_file}[/cyan]")
    
    return r1_file, r2_file, wl_path

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic 10x single-cell FASTQ data")
    parser.add_argument("--outdir", default="data/raw", help="Output directory for FASTQ files")
    parser.add_argument("--whitelist", default="data/whitelist/737K-august-2016.txt", help="Output path for whitelist")
    parser.add_argument("--sample-id", default="scRNA_sample_01", help="Sample identifier prefix")
    parser.add_argument("--num-cells", type=int, default=100, help="Number of real single cells to simulate")
    parser.add_argument("--reads-per-cell", type=int, default=200, help="Average reads per cell")
    parser.add_argument("--num-ambient", type=int, default=300, help="Number of ambient background droplets")
    args = parser.parse_args()
    
    generate_single_cell_dataset(
        outdir=args.outdir,
        whitelist_file=args.whitelist,
        sample_id=args.sample_id,
        num_cells=args.num_cells,
        reads_per_cell=args.reads_per_cell,
        num_ambient_cells=args.num_ambient
    )

if __name__ == "__main__":
    main()
