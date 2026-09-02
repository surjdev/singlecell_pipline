"""
Smart-seq2 & Full-Length scRNA-seq FASTQ Preprocessing Engine.
Performs paired-end trimming via fastp: Nextera adapter trimming,
ISPCR oligo removal, poly-G/poly-X tail trimming, sliding-window quality filtering.
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Tuple
from rich.console import Console
from rich.table import Table

from src.utils import run_cmd, ensure_dir, load_config

console = Console()

class SmartSeq2Preprocessor:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.prep_cfg = config["preprocessing"]
        self.fastp_cfg = self.prep_cfg.get("fastp", {})
        self.threads = self.prep_cfg.get("threads", 4)
        
        self.raw_dir = Path(self.paths["raw_dir"])
        self.clean_dir = ensure_dir(self.paths["clean_dir"])
        self.reports_dir = ensure_dir(self.paths.get("fastp_dir", "reports/fastp"))
        
    def find_cell_samples(self) -> List[Tuple[str, Path, Path]]:
        """Discover paired-end FASTQ samples in raw_dir."""
        r1_files = sorted(list(self.raw_dir.glob("*_R1.fastq.gz")) + list(self.raw_dir.glob("*_1.fastq.gz")))
        samples = []
        for r1 in r1_files:
            if "_R1.fastq.gz" in r1.name:
                cell_id = r1.name.replace("_R1.fastq.gz", "")
                r2 = self.raw_dir / f"{cell_id}_R2.fastq.gz"
            else:
                cell_id = r1.name.replace("_1.fastq.gz", "")
                r2 = self.raw_dir / f"{cell_id}_2.fastq.gz"
                
            if r2.exists():
                samples.append((cell_id, r1, r2))
            else:
                console.print(f"[yellow]Warning: Mate R2 not found for {r1.name}[/yellow]")
        return samples

    def trim_sample(self, cell_id: str, r1_in: Path, r2_in: Path) -> Dict[str, Any]:
        """Run fastp paired-end trimming for a single cell."""
        out_r1 = self.clean_dir / f"{cell_id}_val_R1.fastq.gz"
        out_r2 = self.clean_dir / f"{cell_id}_val_R2.fastq.gz"
        html_report = self.reports_dir / f"{cell_id}_fastp.html"
        json_report = self.reports_dir / f"{cell_id}_fastp.json"
        
        cmd = [
            "fastp",
            "-i", str(r1_in),
            "-I", str(r2_in),
            "-o", str(out_r1),
            "-O", str(out_r2),
            "-h", str(html_report),
            "-j", str(json_report),
            "-w", str(self.threads),
            "-q", str(self.fastp_cfg.get("qualified_quality_phred", 20)),
            "-u", str(self.fastp_cfg.get("unqualified_percent_limit", 30)),
            "-n", str(self.fastp_cfg.get("n_base_limit", 3)),
            "-l", str(self.fastp_cfg.get("min_length", 25))
        ]
        
        # Poly-G trimming
        if self.fastp_cfg.get("trim_poly_g", True):
            cmd.append("--trim_poly_g")
            cmd.extend(["--poly_g_min_len", str(self.fastp_cfg.get("poly_g_min_len", 10))])
            
        # Poly-X trimming
        if self.fastp_cfg.get("trim_poly_x", True):
            cmd.append("--trim_poly_x")
            cmd.extend(["--poly_x_min_len", str(self.fastp_cfg.get("poly_x_min_len", 10))])
            
        # Cut front / tail
        if self.fastp_cfg.get("cut_front", True):
            cmd.append("--cut_front")
            cmd.extend(["--cut_front_window_size", str(self.fastp_cfg.get("cut_front_window_size", 4))])
            cmd.extend(["--cut_front_mean_quality", str(self.fastp_cfg.get("cut_front_mean_quality", 20))])
            
        if self.fastp_cfg.get("cut_tail", True):
            cmd.append("--cut_tail")
            cmd.extend(["--cut_tail_window_size", str(self.fastp_cfg.get("cut_tail_window_size", 4))])
            cmd.extend(["--cut_tail_mean_quality", str(self.fastp_cfg.get("cut_tail_mean_quality", 20))])
            
        # Adapter sequence
        adapter1 = self.fastp_cfg.get("adapter_sequence")
        if adapter1:
            cmd.extend(["--adapter_sequence", adapter1])
        adapter2 = self.fastp_cfg.get("adapter_sequence_r2")
        if adapter2:
            cmd.extend(["--adapter_sequence_r2", adapter2])
            
        run_cmd(cmd, desc=f"fastp trimming [{cell_id}]")
        
        # Parse JSON summary
        metrics = {"cell_id": cell_id}
        if json_report.exists():
            with open(json_report, "r") as f:
                data = json.load(f)
                summary = data.get("summary", {})
                before = summary.get("before_filtering", {})
                after = summary.get("after_filtering", {})
                metrics.update({
                    "raw_reads": before.get("total_reads", 0),
                    "clean_reads": after.get("total_reads", 0),
                    "raw_q30_rate": before.get("q30_rate", 0.0),
                    "clean_q30_rate": after.get("q30_rate", 0.0),
                    "passed_filter_rate": (after.get("total_reads", 0) / max(1, before.get("total_reads", 1))),
                    "adapter_trimmed_reads": data.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
                })
        return metrics

    def run_all(self) -> List[Dict[str, Any]]:
        """Run fastp trimming across all discovered cells."""
        samples = self.find_cell_samples()
        if not samples:
            console.print(f"[yellow]No sample FASTQs found in {self.raw_dir}[/yellow]")
            return []
            
        console.print(f"[bold cyan]Found {len(samples)} cell FASTQ pair(s). Starting fastp preprocessing...[/bold cyan]")
        results = []
        for cell_id, r1, r2 in samples:
            m = self.trim_sample(cell_id, r1, r2)
            results.append(m)
            
        # Print summary table
        table = Table(title="Smart-seq2 fastp Preprocessing Summary", header_style="bold magenta", border_style="cyan")
        table.add_column("Cell ID", style="bold white")
        table.add_column("Raw Reads", justify="right")
        table.add_column("Clean Reads", justify="right")
        table.add_column("Passed Filter", justify="right", style="green")
        table.add_column("Raw Q30", justify="right")
        table.add_column("Clean Q30", justify="right", style="cyan")
        
        for r in results:
            table.add_row(
                r["cell_id"],
                f"{r.get('raw_reads', 0):,}",
                f"{r.get('clean_reads', 0):,}",
                f"{r.get('passed_filter_rate', 0.0)*100:.1f}%",
                f"{r.get('raw_q30_rate', 0.0)*100:.1f}%",
                f"{r.get('clean_q30_rate', 0.0)*100:.1f}%"
            )
        console.print(table)
        return results

def main():
    config = load_config()
    preprocessor = SmartSeq2Preprocessor(config)
    preprocessor.run_all()

if __name__ == "__main__":
    main()
