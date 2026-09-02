"""
Splice-Aware Sequence Aligner Module (STAR & HISAT2).
Handles reference genome indexing, paired-end alignment per cell,
BAM sorting, indexing via samtools, and alignment QC reporting.
"""

import re
import os
import shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
from rich.console import Console
from rich.table import Table

from src.utils import run_cmd, ensure_dir, load_config

console = Console()

class STARAligner:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.align_cfg = config.get("alignment", {})
        self.star_cfg = self.align_cfg.get("star", {})
        self.threads = self.align_cfg.get("threads", 4)
        
        self.clean_dir = Path(self.paths["clean_dir"])
        self.aligned_dir = ensure_dir(self.paths.get("aligned_dir", "data/aligned"))
        self.index_dir = ensure_dir(self.paths["star_index_dir"])
        self.genome_fasta = Path(self.paths["genome_fasta"])
        self.gtf_file = Path(self.paths["gtf_file"])
        
    def build_index_if_needed(self):
        """Generate STAR genome index if SAindex or Genome is missing."""
        sa_file = self.index_dir / "SA"
        if sa_file.exists() and (self.index_dir / "Genome").exists():
            console.print(f"[dim]STAR index already exists at {self.index_dir}[/dim]")
            return
            
        console.print(f"[bold cyan]Building STAR genome index in {self.index_dir}...[/bold cyan]")
        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeDir", str(self.index_dir),
            "--genomeFastaFiles", str(self.genome_fasta),
            "--sjdbGTFfile", str(self.gtf_file),
            "--runThreadN", str(self.threads),
            "--genomeSAindexNbases", str(self.star_cfg.get("genomeSAindexNbases", 6))
        ]
        run_cmd(cmd, desc="STAR genomeGenerate")
        console.print("[green]✔ STAR genome index built successfully![/green]")

    def align_cell(self, cell_id: str, r1_path: Path, r2_path: Path) -> Dict[str, Any]:
        """Align paired-end reads for a single cell using STAR."""
        out_prefix = self.aligned_dir / f"{cell_id}_"
        final_bam = self.aligned_dir / f"{cell_id}_Aligned.sortedByCoord.out.bam"
        log_final = self.aligned_dir / f"{cell_id}_Log.final.out"
        
        cmd = [
            "STAR",
            "--genomeDir", str(self.index_dir),
            "--readFilesIn", str(r1_path), str(r2_path),
            "--readFilesCommand", "zcat",
            "--outFileNamePrefix", str(out_prefix),
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--quantMode", self.star_cfg.get("quantMode", "GeneCounts"),
            "--runThreadN", str(self.threads),
            "--outFilterMultimapNmax", str(self.star_cfg.get("outFilterMultimapNmax", 20)),
            "--outFilterMismatchNmax", str(self.star_cfg.get("outFilterMismatchNmax", 10)),
            "--alignIntronMin", str(self.star_cfg.get("alignIntronMin", 20)),
            "--alignIntronMax", str(self.star_cfg.get("alignIntronMax", 1000000)),
            "--alignMatesGapMax", str(self.star_cfg.get("alignMatesGapMax", 1000000))
        ]
        run_cmd(cmd, desc=f"STAR alignment [{cell_id}]")
        
        # Index BAM with samtools
        if final_bam.exists():
            run_cmd(["samtools", "index", "-@", str(self.threads), str(final_bam)], desc=f"samtools index [{cell_id}]")
            
        # Parse STAR Log.final.out
        stats = {"cell_id": cell_id, "aligner": "STAR"}
        if log_final.exists():
            with open(log_final, "r") as f:
                for line in f:
                    if "Number of input reads" in line:
                        stats["input_reads"] = int(line.split("|")[1].strip())
                    elif "Uniquely mapped reads number" in line:
                        stats["uniquely_mapped"] = int(line.split("|")[1].strip())
                    elif "Uniquely mapped reads %" in line:
                        stats["uniquely_mapped_pct"] = float(line.split("|")[1].strip().replace("%", ""))
                    elif "% of reads mapped to multiple loci" in line:
                        stats["multimapped_pct"] = float(line.split("|")[1].strip().replace("%", ""))
                    elif "% of reads unmapped: too short" in line:
                        stats["unmapped_pct"] = float(line.split("|")[1].strip().replace("%", ""))
        return stats

    def run_all(self) -> List[Dict[str, Any]]:
        self.build_index_if_needed()
        r1_files = sorted(list(self.clean_dir.glob("*_val_R1.fastq.gz")) + list(self.clean_dir.glob("*_val_1.fastq.gz")))
        if not r1_files:
            # Fallback to clean dir any fastqs
            r1_files = sorted(list(self.clean_dir.glob("*_R1.fastq.gz")))
            
        results = []
        for r1 in r1_files:
            cell_id = r1.name.replace("_val_R1.fastq.gz", "").replace("_val_1.fastq.gz", "").replace("_R1.fastq.gz", "")
            r2 = self.clean_dir / r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz").replace("_1.fastq.gz", "_2.fastq.gz")
            if r2.exists():
                st = self.align_cell(cell_id, r1, r2)
                results.append(st)
                
        self._print_summary(results)
        return results

    def _print_summary(self, results: List[Dict[str, Any]]):
        table = Table(title="STAR Alignment Summary", header_style="bold magenta", border_style="cyan")
        table.add_column("Cell ID", style="bold white")
        table.add_column("Input Reads", justify="right")
        table.add_column("Uniquely Mapped", justify="right", style="green")
        table.add_column("Unique %", justify="right", style="bold green")
        table.add_column("Multi %", justify="right", style="yellow")
        
        for r in results:
            table.add_row(
                r["cell_id"],
                f"{r.get('input_reads', 0):,}",
                f"{r.get('uniquely_mapped', 0):,}",
                f"{r.get('uniquely_mapped_pct', 0.0):.1f}%",
                f"{r.get('multimapped_pct', 0.0):.1f}%"
            )
        console.print(table)


class HISAT2Aligner:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.align_cfg = config.get("alignment", {})
        self.threads = self.align_cfg.get("threads", 4)
        
        self.clean_dir = Path(self.paths["clean_dir"])
        self.aligned_dir = ensure_dir(self.paths.get("aligned_dir", "data/aligned"))
        self.index_prefix = self.paths["hisat2_index_prefix"]
        self.genome_fasta = Path(self.paths["genome_fasta"])
        ensure_dir(Path(self.index_prefix).parent)

    def build_index_if_needed(self):
        """Build HISAT2 index if not present."""
        test_file = Path(f"{self.index_prefix}.1.ht2")
        if test_file.exists():
            console.print(f"[dim]HISAT2 index already exists at {self.index_prefix}[/dim]")
            return
            
        console.print(f"[bold cyan]Building HISAT2 index at {self.index_prefix}...[/bold cyan]")
        cmd = [
            "hisat2-build",
            "-p", str(self.threads),
            str(self.genome_fasta),
            str(self.index_prefix)
        ]
        run_cmd(cmd, desc="hisat2-build")
        console.print("[green]✔ HISAT2 index built successfully![/green]")

    def align_cell(self, cell_id: str, r1_path: Path, r2_path: Path) -> Dict[str, Any]:
        """Align paired-end reads using HISAT2 and stream to sorted BAM."""
        out_bam = self.aligned_dir / f"{cell_id}_hisat2_sorted.bam"
        log_file = self.aligned_dir / f"{cell_id}_hisat2.log"
        
        pipe_cmd = (
            f"hisat2 -p {self.threads} --dta -x {self.index_prefix} "
            f"-1 {r1_path} -2 {r2_path} 2> {log_file} | "
            f"samtools sort -@ {self.threads} -o {out_bam} -"
        )
        run_cmd(pipe_cmd, desc=f"HISAT2 alignment [{cell_id}]")
        
        # Index BAM
        if out_bam.exists():
            run_cmd(["samtools", "index", "-@", str(self.threads), str(out_bam)], desc=f"samtools index [{cell_id}]")
            
        # Parse alignment rate from log
        stats = {"cell_id": cell_id, "aligner": "HISAT2"}
        if log_file.exists():
            with open(log_file, "r") as f:
                content = f.read()
                m_rate = re.search(r"([\d\.]+)% overall alignment rate", content)
                if m_rate:
                    stats["alignment_rate"] = float(m_rate.group(1))
        return stats

    def run_all(self) -> List[Dict[str, Any]]:
        self.build_index_if_needed()
        r1_files = sorted(list(self.clean_dir.glob("*_val_R1.fastq.gz")) + list(self.clean_dir.glob("*_val_1.fastq.gz")))
        if not r1_files:
            r1_files = sorted(list(self.clean_dir.glob("*_R1.fastq.gz")))
            
        results = []
        for r1 in r1_files:
            cell_id = r1.name.replace("_val_R1.fastq.gz", "").replace("_val_1.fastq.gz", "").replace("_R1.fastq.gz", "")
            r2 = self.clean_dir / r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz").replace("_1.fastq.gz", "_2.fastq.gz")
            if r2.exists():
                st = self.align_cell(cell_id, r1, r2)
                results.append(st)
        return results

def main():
    config = load_config()
    aligner_name = config.get("alignment", {}).get("default_aligner", "star").lower()
    if aligner_name == "hisat2":
        aligner = HISAT2Aligner(config)
    else:
        aligner = STARAligner(config)
    aligner.run_all()

if __name__ == "__main__":
    main()
