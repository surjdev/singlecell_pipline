"""
Smart-seq2 & Transcriptomic Quality Control Engine.
Handles FastQC execution on multi-sample paired FASTQs,
MultiQC dashboard compilation, and quality metric extraction.
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
from rich.console import Console

from src.utils import run_cmd, ensure_dir, load_config

console = Console()

class SmartSeq2QC:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.threads = config.get("preprocessing", {}).get("threads", 4)
        
        self.raw_dir = Path(self.paths["raw_dir"])
        self.clean_dir = Path(self.paths["clean_dir"])
        self.qc_raw_dir = ensure_dir(self.paths.get("qc_raw_dir", "reports/qc_raw"))
        self.qc_clean_dir = ensure_dir(self.paths.get("qc_clean_dir", "reports/qc_clean"))
        self.multiqc_dir = ensure_dir(self.paths.get("multiqc_dir", "reports/multiqc"))
        self.reports_dir = ensure_dir(self.paths.get("reports_dir", "reports"))
        
    def run_fastqc_stage(self, stage: str = "raw") -> Path:
        """Run FastQC across all FASTQ files in raw or clean directory."""
        if stage == "raw":
            in_dir = self.raw_dir
            out_dir = self.qc_raw_dir
        else:
            in_dir = self.clean_dir
            out_dir = self.qc_clean_dir
            
        fastq_files = sorted(list(in_dir.glob("*.fastq.gz")) + list(in_dir.glob("*.fq.gz")))
        if not fastq_files:
            console.print(f"[yellow]No FASTQ files found in {in_dir} for FastQC.[/yellow]")
            return out_dir
            
        console.print(f"[bold cyan]Running FastQC on {len(fastq_files)} file(s) for stage: {stage}...[/bold cyan]")
        cmd = [
            "fastqc",
            "-o", str(out_dir),
            "-t", str(self.threads),
            "--quiet"
        ] + [str(f) for f in fastq_files]
        
        run_cmd(cmd, desc=f"FastQC ({stage})")
        console.print(f"[green]✔ FastQC complete. Reports saved to {out_dir}[/green]")
        return out_dir

    def run_multiqc(self) -> Path:
        """Compile MultiQC report aggregating FastQC, fastp, STAR/HISAT2, and featureCounts."""
        console.print(f"[bold cyan]Aggregating all QC logs into MultiQC dashboard...[/bold cyan]")
        
        # Scan reports_dir, data/aligned, data/counts
        search_dirs = [
            str(self.reports_dir),
            str(self.paths.get("aligned_dir", "data/aligned")),
            str(self.paths.get("counts_dir", "data/counts"))
        ]
        
        cmd = [
            "multiqc",
            "--outdir", str(self.multiqc_dir),
            "--force",
            "--title", self.config.get("project", {}).get("name", "Smart-seq2 Pipeline"),
            "--filename", "multiqc_report.html"
        ] + [d for d in search_dirs if Path(d).exists()]
        
        run_cmd(cmd, desc="MultiQC report compilation")
        html_report = self.multiqc_dir / "multiqc_report.html"
        if html_report.exists():
            console.print(f"[bold green]✔ MultiQC Dashboard compiled successfully:[/bold green] {html_report}")
        return self.multiqc_dir

def main():
    config = load_config()
    qc = SmartSeq2QC(config)
    qc.run_fastqc_stage("raw")
    qc.run_multiqc()

if __name__ == "__main__":
    main()
