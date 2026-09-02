"""
Synthetic Data Generator for Smart-seq2 Single-Cell RNA-seq Pipeline.
Generates a mini reference genome, GTF annotation file, and realistic
per-cell paired-end FASTQ files with Nextera adapters, sequencing errors,
and variable transcript expression across cells.
"""

import os
import gzip
import random
import argparse
from pathlib import Path
from typing import List, Dict, Tuple
from rich.console import Console

from src.utils import ensure_dir

console = Console()

# Fixed random seed for reproducibility
random.seed(42)

BASES = ["A", "C", "G", "T"]
NEXTERA_ADAPTER = "CTGTCTCTTATACACATCT"
TSO_SEQ = "AAGCAGTGGTATCAACGCAGAGTACATGGG"

def generate_random_seq(length: int) -> str:
    return "".join(random.choices(BASES, k=length))

def reverse_complement(seq: str) -> str:
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]

def make_phred_string(length: int, base_q: int = 37, decay: bool = True) -> str:
    """Generate Phred quality scores with optional 3' quality decay."""
    scores = []
    for i in range(length):
        if decay:
            # Gradual drop towards the end with random noise
            q = max(10, int(base_q - (i / length) * 8 + random.randint(-4, 3)))
        else:
            q = max(20, int(base_q + random.randint(-3, 2)))
        q = min(40, max(2, q))
        scores.append(chr(q + 33))
    return "".join(scores)

def create_mini_reference(ref_dir: Path) -> Tuple[Path, Path, List[Dict]]:
    """
    Creates a synthetic mini genome FASTA and corresponding GTF file
    containing multiple genes with multi-exon structures.
    """
    ensure_dir(ref_dir)
    fasta_path = ref_dir / "mini_genome.fa"
    gtf_path = ref_dir / "mini_genes.gtf"
    
    genes = []
    chrom_seq_parts = []
    current_pos = 1000
    
    # 15 synthetic genes (including a couple housekeeping & variable genes)
    gene_names = [
        "GAPDH_like", "ACTB_like", "TP53_like", "CD4_like", "CD8A_like",
        "MS4A1_like", "GNLY_like", "NKG7_like", "IL7R_like", "MT-ND1_like",
        "MT-CO1_like", "MALAT1_like", "SOX2_like", "OCT4_like", "NANOG_like"
    ]
    
    gtf_lines = [
        "##gff-version 2",
        "#!genome-build mini_genome_v1.0",
        "#!genome-version mini_genome_v1.0"
    ]
    
    for i, gname in enumerate(gene_names, 1):
        gid = f"GENE_{i:03d}_{gname}"
        tid = f"TRANS_{i:03d}_1"
        strand = "+" if i % 2 == 1 else "-"
        
        # 2 to 4 exons per gene
        num_exons = random.randint(2, 4)
        exon_lengths = [random.randint(150, 350) for _ in range(num_exons)]
        intron_lengths = [random.randint(100, 250) for _ in range(num_exons - 1)]
        
        gene_start = current_pos
        exons_info = []
        transcript_seq_parts = []
        
        exon_start = gene_start
        for e_idx, e_len in enumerate(exon_lengths):
            exon_end = exon_start + e_len - 1
            e_seq = generate_random_seq(e_len)
            transcript_seq_parts.append(e_seq)
            chrom_seq_parts.append(e_seq)
            
            exons_info.append((exon_start, exon_end, e_idx + 1))
            
            if e_idx < len(intron_lengths):
                # Add intron sequence
                i_len = intron_lengths[e_idx]
                chrom_seq_parts.append(generate_random_seq(i_len))
                exon_start = exon_end + 1 + i_len
            else:
                exon_start = exon_end + 1
                
        gene_end = exon_end
        current_pos = gene_end + random.randint(300, 600)
        
        full_tx_seq = "".join(transcript_seq_parts)
        if strand == "-":
            full_tx_seq = reverse_complement(full_tx_seq)
            
        genes.append({
            "gene_id": gid,
            "gene_name": gname,
            "transcript_id": tid,
            "strand": strand,
            "start": gene_start,
            "end": gene_end,
            "seq": full_tx_seq,
            "exons": exons_info
        })
        
        # GTF gene entry
        gtf_lines.append(
            f'chr1\tSynthetic\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t'
            f'gene_id "{gid}"; gene_name "{gname}";'
        )
        # GTF transcript entry
        gtf_lines.append(
            f'chr1\tSynthetic\ttranscript\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t'
            f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}";'
        )
        # GTF exon entries
        for estart, eend, enum in exons_info:
            gtf_lines.append(
                f'chr1\tSynthetic\texon\t{estart}\t{eend}\t.\t{strand}\t.\t'
                f'gene_id "{gid}"; transcript_id "{tid}"; exon_number "{enum}"; gene_name "{gname}";'
            )
            
    # Assemble chromosome sequence
    chrom1_seq = "".join(chrom_seq_parts)
    
    # Write FASTA
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        # 80 chars per line
        for j in range(0, len(chrom1_seq), 80):
            f.write(chrom1_seq[j:j+80] + "\n")
            
    # Write GTF
    with open(gtf_path, "w") as f:
        f.write("\n".join(gtf_lines) + "\n")
        
    console.print(f"[green]✔ Generated mini reference genome:[/green] {fasta_path} ({len(chrom1_seq):,} bp)")
    console.print(f"[green]✔ Generated mini GTF annotation:[/green] {gtf_path} ({len(genes)} genes)")
    
    return fasta_path, gtf_path, genes

def generate_smartseq2_cell_fastqs(
    genes: List[Dict],
    out_dir: Path,
    num_cells: int = 8,
    reads_per_cell: int = 500,
    read_len: int = 100
) -> List[Tuple[str, Path, Path]]:
    """
    Generates paired-end FASTQs for each individual cell (e.g. cell_01_R1.fastq.gz, cell_01_R2.fastq.gz).
    Simulates:
    - Differing cell expression profiles (cell states / types)
    - Fragment size distribution (insert size ~250-400bp)
    - Nextera adapter read-through contamination on short inserts (~5% reads)
    - Sequencing errors / low-quality 3' tails
    """
    ensure_dir(out_dir)
    samples = []
    
    for cell_idx in range(1, num_cells + 1):
        cell_id = f"cell_{cell_idx:02d}"
        r1_path = out_dir / f"{cell_id}_R1.fastq.gz"
        r2_path = out_dir / f"{cell_id}_R2.fastq.gz"
        
        # Define expression weights for this cell
        weights = []
        for g in genes:
            gname = g["gene_name"]
            # Housekeeping genes express consistently
            if "GAPDH" in gname or "ACTB" in gname:
                w = random.uniform(8.0, 15.0)
            elif "CD4" in gname or "IL7R" in gname:
                w = random.uniform(5.0, 12.0) if cell_idx <= num_cells // 2 else random.uniform(0.1, 1.0)
            elif "CD8A" in gname or "NKG7" in gname:
                w = random.uniform(5.0, 12.0) if cell_idx > num_cells // 2 else random.uniform(0.1, 1.0)
            else:
                w = random.uniform(0.5, 6.0)
            weights.append(w)
            
        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for r_i in range(reads_per_cell):
                chosen_gene = random.choices(genes, weights=weights, k=1)[0]
                tx_seq = chosen_gene["seq"]
                tx_len = len(tx_seq)
                
                # Pick fragment
                insert_size = random.randint(150, min(500, tx_len))
                if tx_len <= insert_size:
                    frag_seq = tx_seq
                else:
                    start = random.randint(0, tx_len - insert_size)
                    frag_seq = tx_seq[start:start + insert_size]
                    
                frag_len = len(frag_seq)
                
                # R1 starts from 5' of fragment
                if frag_len < read_len:
                    # Adapter read-through
                    r1_seq = frag_seq + NEXTERA_ADAPTER[:read_len - frag_len]
                else:
                    r1_seq = frag_seq[:read_len]
                    
                # R2 starts from 3' of fragment (reverse complement)
                rc_frag = reverse_complement(frag_seq)
                if frag_len < read_len:
                    r2_seq = rc_frag + NEXTERA_ADAPTER[:read_len - frag_len]
                else:
                    r2_seq = rc_frag[:read_len]
                    
                # Random sequencing error (~0.5%)
                r1_list = list(r1_seq)
                for pos in range(len(r1_list)):
                    if random.random() < 0.005:
                        r1_list[pos] = random.choice(BASES)
                r1_seq = "".join(r1_list)
                
                r2_list = list(r2_seq)
                for pos in range(len(r2_list)):
                    if random.random() < 0.005:
                        r2_list[pos] = random.choice(BASES)
                r2_seq = "".join(r2_list)
                
                # Write FASTQ records
                read_name = f"@{cell_id}_READ_{r_i+1}"
                f1.write(f"{read_name}/1\n{r1_seq}\n+\n{make_phred_string(len(r1_seq))}\n")
                f2.write(f"{read_name}/2\n{r2_seq}\n+\n{make_phred_string(len(r2_seq))}\n")
                
        samples.append((cell_id, r1_path, r2_path))
        console.print(f"[green]✔ Generated Smart-seq2 FASTQ pair:[/green] {cell_id} ({reads_per_cell:,} read pairs)")
        
    return samples

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic Smart-seq2 FASTQ data & mini reference genome")
    parser.add_argument("--num-cells", type=int, default=8, help="Number of cells/wells to generate")
    parser.add_argument("--reads-per-cell", type=int, default=500, help="Number of read pairs per cell")
    parser.add_argument("--ref-dir", type=str, default="data/reference", help="Reference output directory")
    parser.add_argument("--raw-dir", type=str, default="data/raw", help="Raw FASTQ output directory")
    
    args = parser.parse_args()
    
    ref_dir = Path(args.ref_dir)
    raw_dir = Path(args.raw_dir)
    
    console.print("[bold blue]Generating Smart-seq2 Synthetic Test Dataset...[/bold blue]")
    _, _, genes = create_mini_reference(ref_dir)
    generate_smartseq2_cell_fastqs(
        genes=genes,
        out_dir=raw_dir,
        num_cells=args.num_cells,
        reads_per_cell=args.reads_per_cell
    )
    console.print("[bold green]✔ Synthetic Smart-seq2 dataset and reference created successfully![/bold green]")

if __name__ == "__main__":
    main()
