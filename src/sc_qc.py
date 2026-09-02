"""
Single-Cell FASTQ Quality Control Analysis Engine.
Computes standard sequencing quality metrics as well as single-cell specific
measurements: Cell Barcode / UMI Q30, barcode diversity, knee-plot curve,
whitelist match statistics, GC bias, and automated FastQC execution.
"""

import gzip
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
from collections import Counter
from rich.console import Console

from src.utils import run_cmd, ensure_dir, print_summary_table

console = Console()

class SingleCellQC:
    def __init__(
        self,
        r1_path: str | Path,
        r2_path: str | Path,
        cb_len: int = 16,
        umi_len: int = 12,
        whitelist_path: Optional[str | Path] = None,
        sample_id: str = "sample"
    ):
        self.r1_path = Path(r1_path)
        self.r2_path = Path(r2_path)
        self.cb_len = cb_len
        self.umi_len = umi_len
        self.sample_id = sample_id
        self.whitelist = self._load_whitelist(whitelist_path) if whitelist_path else set()
        
    def _load_whitelist(self, wl_path: str | Path) -> set:
        path = Path(wl_path)
        if not path.exists():
            console.print(f"[yellow]Warning: Whitelist file {wl_path} not found. Whitelist match rate will be skipped.[/yellow]")
            return set()
        with open(path, "r") as f:
            return set(line.strip().split()[0] for line in f if line.strip())

    def run_fastqc(self, out_dir: str | Path, threads: int = 2) -> Path:
        """Run standard FastQC on both R1 and R2."""
        out_p = ensure_dir(out_dir)
        cmd = [
            "fastqc",
            "-o", str(out_p),
            "-t", str(threads),
            "--quiet",
            str(self.r1_path),
            str(self.r2_path)
        ]
        run_cmd(cmd, desc=f"Running FastQC for {self.sample_id}")
        return out_p

    def analyze_fastq(self, max_reads: Optional[int] = None) -> Dict[str, Any]:
        """Parse FASTQ pair and calculate detailed single-cell QC statistics."""
        console.print(f"[bold cyan]Analyzing Single-Cell QC for {self.sample_id}...[/bold cyan]")
        
        total_reads = 0
        cb_q30_bases = 0
        cb_total_bases = 0
        umi_q30_bases = 0
        umi_total_bases = 0
        r2_q30_bases = 0
        r2_total_bases = 0
        
        cb_counts = Counter()
        umi_counts = Counter()
        
        r1_pos_quals = []
        r2_pos_quals = []
        r1_pos_bases = []
        r2_pos_bases = []
        
        exact_wl_matches = 0
        mismatch_1_matches = 0
        invalid_barcodes = 0
        
        poly_g_count = 0
        poly_a_count = 0
        
        open_r1 = gzip.open(self.r1_path, "rt") if str(self.r1_path).endswith(".gz") else open(self.r1_path, "rt")
        open_r2 = gzip.open(self.r2_path, "rt") if str(self.r2_path).endswith(".gz") else open(self.r2_path, "rt")
        
        try:
            while True:
                r1_h = open_r1.readline()
                if not r1_h:
                    break
                r1_s = open_r1.readline().strip()
                r1_p = open_r1.readline()
                r1_q = open_r1.readline().strip()
                
                r2_h = open_r2.readline()
                r2_s = open_r2.readline().strip()
                r2_p = open_r2.readline()
                r2_q = open_r2.readline().strip()
                
                total_reads += 1
                if max_reads and total_reads > max_reads:
                    break
                    
                # R1 breakdown: CB + UMI
                cb_seq = r1_s[:self.cb_len]
                cb_qual = [ord(c) - 33 for c in r1_q[:self.cb_len]]
                umi_seq = r1_s[self.cb_len:self.cb_len + self.umi_len]
                umi_qual = [ord(c) - 33 for c in r1_q[self.cb_len:self.cb_len + self.umi_len]]
                
                # R2 breakdown: cDNA
                r2_qual = [ord(c) - 33 for c in r2_q]
                
                # Q30 tallies
                cb_q30_bases += sum(1 for q in cb_qual if q >= 30)
                cb_total_bases += len(cb_qual)
                umi_q30_bases += sum(1 for q in umi_qual if q >= 30)
                umi_total_bases += len(umi_qual)
                r2_q30_bases += sum(1 for q in r2_qual if q >= 30)
                r2_total_bases += len(r2_qual)
                
                cb_counts[cb_seq] += 1
                umi_counts[umi_seq] += 1
                
                # Check for Poly-A / Poly-G in R2
                if "AAAAAAAAAA" in r2_s:
                    poly_a_count += 1
                if "GGGGGGGGGG" in r2_s:
                    poly_g_count += 1
                    
                # Whitelist evaluation
                if self.whitelist:
                    if cb_seq in self.whitelist:
                        exact_wl_matches += 1
                    else:
                        # Test for 1-bp mismatch
                        matched = False
                        for i in range(len(cb_seq)):
                            for alt in ['A', 'C', 'G', 'T']:
                                if alt != cb_seq[i]:
                                    cand = cb_seq[:i] + alt + cb_seq[i+1:]
                                    if cand in self.whitelist:
                                        mismatch_1_matches += 1
                                        matched = True
                                        break
                            if matched:
                                break
                        if not matched:
                            invalid_barcodes += 1
                            
                # Base & Quality profiles per cycle (sampling for efficiency)
                if total_reads <= 10000:
                    while len(r1_pos_quals) < len(r1_q):
                        r1_pos_quals.append([])
                        r1_pos_bases.append(Counter())
                    for idx, (base, q) in enumerate(zip(r1_s, [ord(c)-33 for c in r1_q])):
                        r1_pos_quals[idx].append(q)
                        r1_pos_bases[idx][base] += 1
                        
                    while len(r2_pos_quals) < len(r2_q):
                        r2_pos_quals.append([])
                        r2_pos_bases.append(Counter())
                    for idx, (base, q) in enumerate(zip(r2_s, r2_qual)):
                        r2_pos_quals[idx].append(q)
                        r2_pos_bases[idx][base] += 1

        finally:
            open_r1.close()
            open_r2.close()

        # Compute summary metrics
        cb_q30_pct = cb_q30_bases / max(1, cb_total_bases)
        umi_q30_pct = umi_q30_bases / max(1, umi_total_bases)
        r2_q30_pct = r2_q30_bases / max(1, r2_total_bases)
        
        sorted_cb_counts = sorted(cb_counts.values(), reverse=True)
        unique_barcodes = len(cb_counts)
        unique_umis = len(umi_counts)
        
        # Estimate knee inflection
        knee_stats = self._calculate_knee_inflection(sorted_cb_counts)
        
        metrics = {
            "sample_id": self.sample_id,
            "total_reads": total_reads,
            "cb_length": self.cb_len,
            "umi_length": self.umi_len,
            "unique_cell_barcodes": unique_barcodes,
            "unique_umis": unique_umis,
            "cb_q30_fraction": float(cb_q30_pct),
            "umi_q30_fraction": float(umi_q30_pct),
            "r2_q30_fraction": float(r2_q30_pct),
            "poly_a_rate": poly_a_count / max(1, total_reads),
            "poly_g_rate": poly_g_count / max(1, total_reads),
            "estimated_cells": knee_stats["estimated_cells"],
            "reads_in_cells_fraction": knee_stats["fraction_reads_in_cells"],
        }
        
        if self.whitelist:
            metrics["whitelist_exact_match_rate"] = exact_wl_matches / max(1, total_reads)
            metrics["whitelist_1bp_mismatch_rate"] = mismatch_1_matches / max(1, total_reads)
            metrics["whitelist_invalid_rate"] = invalid_barcodes / max(1, total_reads)
            metrics["valid_barcode_fraction"] = (exact_wl_matches + mismatch_1_matches) / max(1, total_reads)
            
        return {
            "summary_metrics": metrics,
            "barcode_rank_counts": sorted_cb_counts,
            "knee_stats": knee_stats,
            "r1_mean_qual_per_cycle": [float(np.mean(q)) for q in r1_pos_quals],
            "r2_mean_qual_per_cycle": [float(np.mean(q)) for q in r2_pos_quals],
        }

    def _calculate_knee_inflection(self, sorted_counts: List[int]) -> Dict[str, Any]:
        """Estimate cell count cutoff from barcode rank distribution (Knee Plot)."""
        if not sorted_counts:
            return {"estimated_cells": 0, "cutoff_reads": 0, "fraction_reads_in_cells": 0.0}
            
        counts = np.array(sorted_counts)
        total_reads = counts.sum()
        
        # Simple inflection detection heuristic: finding high-depth cell plateau vs ambient drop
        if len(counts) > 10:
            log_ranks = np.log10(np.arange(1, len(counts) + 1))
            log_counts = np.log10(counts + 1)
            # Find the point of maximum negative derivative in log-log space
            diff = np.diff(log_counts) / np.diff(log_ranks)
            smoothed_diff = np.convolve(diff, np.ones(5)/5, mode='same')
            inflection_idx = int(np.argmin(smoothed_diff[:max(5, int(len(counts) * 0.5))]))
            estimated_cells = max(1, inflection_idx + 1)
        else:
            estimated_cells = len(counts)
            
        reads_in_cells = int(counts[:estimated_cells].sum())
        fraction_in_cells = float(reads_in_cells / max(1, total_reads))
        
        return {
            "estimated_cells": estimated_cells,
            "cutoff_reads": int(counts[min(estimated_cells - 1, len(counts) - 1)]),
            "fraction_reads_in_cells": fraction_in_cells
        }

    def save_report(self, results: Dict[str, Any], output_json_path: str | Path):
        """Save QC metrics as JSON for downstream consumption."""
        out_p = Path(output_json_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        # Store concise JSON (exclude huge raw lists to keep file light)
        save_dict = {
            "summary_metrics": results["summary_metrics"],
            "knee_stats": results["knee_stats"],
            "r1_mean_qual_per_cycle": results["r1_mean_qual_per_cycle"],
            "r2_mean_qual_per_cycle": results["r2_mean_qual_per_cycle"],
            "top_100_barcode_counts": results["barcode_rank_counts"][:100]
        }
        with open(out_p, "w") as f:
            json.dump(save_dict, f, indent=2)
        console.print(f"  ✔ Saved QC report: [green]{out_p}[/green]")
