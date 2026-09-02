"""
Single-Cell FASTQ Preprocessing Engine.
Performs:
1. Fastp adapter trimming, poly-G / poly-A removal, and quality filtering.
2. Barcode error-correction against 10x whitelist (1-Hamming distance).
3. Header extraction or paired cleaned FASTQ generation.
"""

import gzip
import json
import shutil
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple
from collections import defaultdict
from rich.console import Console

from src.utils import run_cmd, ensure_dir, print_summary_table

console = Console()

class SingleCellPreprocessor:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.prep_cfg = config.get("preprocessing", {})
        self.chem_cfg = config.get("chemistry", {})
        
        self.cb_len = self.chem_cfg.get("r1_structure", {}).get("cell_barcode_len", 16)
        self.umi_len = self.chem_cfg.get("r1_structure", {}).get("umi_len", 12)
        
    def _build_whitelist_lookup(self, whitelist_file: str | Path) -> Tuple[Set[str], Dict[str, str]]:
        """
        Build exact whitelist set and 1-Hamming distance error correction lookup.
        """
        exact_wl = set()
        correction_map = {}
        conflict_set = set()
        
        wl_p = Path(whitelist_file)
        if not wl_p.exists():
            console.print(f"[yellow]Warning: Whitelist file {whitelist_file} not found. Error correction disabled.[/yellow]")
            return exact_wl, correction_map
            
        with open(wl_p, "r") as f:
            for line in f:
                bc = line.strip().split()[0]
                if bc:
                    exact_wl.add(bc)
                    
        console.print(f"Building 1-Hamming distance index for {len(exact_wl):,} barcodes...")
        # Precompute 1-bp mutations for fast O(1) correction
        for true_bc in exact_wl:
            for i in range(len(true_bc)):
                for alt in ['A', 'C', 'G', 'T']:
                    if alt != true_bc[i]:
                        mutated = true_bc[:i] + alt + true_bc[i+1:]
                        if mutated in exact_wl:
                            continue # True barcode collision
                        if mutated in correction_map:
                            # Ambiguous 1-bp neighbor pointing to two true barcodes -> mark as conflict
                            conflict_set.add(mutated)
                        else:
                            correction_map[mutated] = true_bc
                            
        # Remove ambiguous conflicts
        for conf in conflict_set:
            if conf in correction_map:
                del correction_map[conf]
                
        console.print(f"  ✔ Whitelist loaded: {len(exact_wl):,} exact barcodes, {len(correction_map):,} correctable 1-bp mutants indexed.")
        return exact_wl, correction_map

    def run_fastp_trimming(
        self,
        r1_in: str | Path,
        r2_in: str | Path,
        r1_out: str | Path,
        r2_out: str | Path,
        report_html: str | Path,
        report_json: str | Path
    ) -> Dict[str, Any]:
        """
        Run fastp for adapter trimming, poly-G/poly-X removal, and quality filtering.
        """
        threads = self.prep_cfg.get("threads", 4)
        fp_cfg = self.prep_cfg.get("fastp", {})
        
        cmd = [
            "fastp",
            "-i", str(r1_in),
            "-I", str(r2_in),
            "-o", str(r1_out),
            "-O", str(r2_out),
            "--html", str(report_html),
            "--json", str(report_json),
            "--thread", str(threads),
            f"--qualified_quality_phred={fp_cfg.get('qualified_quality_phred', 20)}",
            f"--unqualified_percent_limit={fp_cfg.get('unqualified_percent_limit', 30)}",
            f"--n_base_limit={fp_cfg.get('n_base_limit', 3)}",
            f"--length_required={fp_cfg.get('min_length', 25)}",
        ]
        
        if fp_cfg.get("trim_poly_g", True):
            cmd.append("--trim_poly_g")
            cmd.append(f"--poly_g_min_len={fp_cfg.get('poly_g_min_len', 10)}")
            
        if fp_cfg.get("trim_poly_x", True):
            cmd.append("--trim_poly_x")
            cmd.append(f"--poly_x_min_len={fp_cfg.get('poly_x_min_len', 10)}")
            
        if fp_cfg.get("cut_front", True):
            cmd.append("--cut_front")
            cmd.append(f"--cut_front_window_size={fp_cfg.get('cut_front_window_size', 4)}")
            cmd.append(f"--cut_front_mean_quality={fp_cfg.get('cut_front_mean_quality', 20)}")
            
        if fp_cfg.get("cut_tail", True):
            cmd.append("--cut_tail")
            cmd.append(f"--cut_tail_window_size={fp_cfg.get('cut_tail_window_size', 4)}")
            cmd.append(f"--cut_tail_mean_quality={fp_cfg.get('cut_tail_mean_quality', 20)}")
            
        adapter = fp_cfg.get("adapter_sequence_r2")
        if adapter:
            cmd.append(f"--adapter_sequence_r2={adapter}")
        else:
            cmd.append("--disable_adapter_trimming")
        
        ensure_dir(Path(r1_out).parent)
        ensure_dir(Path(report_html).parent)
        
        run_cmd(cmd, desc="Running Fastp Adapter & Quality Trimming")
        
        # Parse fastp JSON report
        with open(report_json, "r") as f:
            fp_stats = json.load(f)
            
        return fp_stats

    def filter_and_correct_barcodes(
        self,
        r1_in: str | Path,
        r2_in: str | Path,
        r1_out: str | Path,
        r2_out: str | Path,
        extracted_out: Optional[str | Path] = None,
        whitelist_file: Optional[str | Path] = None
    ) -> Dict[str, Any]:
        """
        Process trimmed paired FASTQ:
        1. Validate CB against whitelist
        2. Correct 1-bp mismatch CB
        3. Write filtered clean paired FASTQ and/or extracted single FASTQ (cDNA with @NAME_CB_UMI header)
        """
        wl_path = whitelist_file or self.paths.get("whitelist_file")
        exact_wl, correction_map = self._build_whitelist_lookup(wl_path)
        
        console.print(f"[bold cyan]Filtering and error-correcting cell barcodes...[/bold cyan]")
        
        ensure_dir(Path(r1_out).parent)
        
        open_r1_in = gzip.open(r1_in, "rt") if str(r1_in).endswith(".gz") else open(r1_in, "rt")
        open_r2_in = gzip.open(r2_in, "rt") if str(r2_in).endswith(".gz") else open(r2_in, "rt")
        
        open_r1_out = gzip.open(r1_out, "wt")
        open_r2_out = gzip.open(r2_out, "wt")
        open_ext_out = gzip.open(extracted_out, "wt") if extracted_out else None
        
        total_reads = 0
        passed_exact = 0
        passed_corrected = 0
        discarded_invalid = 0
        discarded_short = 0
        
        try:
            while True:
                r1_h = open_r1_in.readline()
                if not r1_h:
                    break
                r1_s = open_r1_in.readline().strip()
                r1_p = open_r1_in.readline()
                r1_q = open_r1_in.readline().strip()
                
                r2_h = open_r2_in.readline()
                r2_s = open_r2_in.readline().strip()
                r2_p = open_r2_in.readline()
                r2_q = open_r2_in.readline().strip()
                
                total_reads += 1
                
                # Check R1 length for full CB + UMI
                req_r1_len = self.cb_len + self.umi_len
                if len(r1_s) < req_r1_len:
                    discarded_short += 1
                    continue
                    
                cb = r1_s[:self.cb_len]
                umi = r1_s[self.cb_len:req_r1_len]
                cb_q = r1_q[:self.cb_len]
                umi_q = r1_q[self.cb_len:req_r1_len]
                
                corrected_cb = cb
                is_valid = True
                
                if exact_wl:
                    if cb in exact_wl:
                        passed_exact += 1
                    elif cb in correction_map:
                        corrected_cb = correction_map[cb]
                        passed_corrected += 1
                    else:
                        discarded_invalid += 1
                        is_valid = False
                else:
                    passed_exact += 1
                    
                if not is_valid:
                    continue
                    
                # Write to clean paired FASTQ
                new_r1_seq = corrected_cb + umi
                new_r1_qual = cb_q + umi_q
                
                open_r1_out.write(f"{r1_h.strip()}\n{new_r1_seq}\n+\n{new_r1_qual}\n")
                open_r2_out.write(f"{r2_h.strip()}\n{r2_s}\n+\n{r2_q}\n")
                
                # Write to extracted FASTQ if requested (standard UMI-tools / Kallisto format)
                if open_ext_out:
                    read_id = r2_h.strip().split()[0]
                    ext_header = f"{read_id}_{corrected_cb}_{umi} {r2_h.strip().split()[1] if len(r2_h.strip().split()) > 1 else '1:N:0:0'}"
                    open_ext_out.write(f"{ext_header}\n{r2_s}\n+\n{r2_q}\n")

        finally:
            open_r1_in.close()
            open_r2_in.close()
            open_r1_out.close()
            open_r2_out.close()
            if open_ext_out:
                open_ext_out.close()
                
        total_kept = passed_exact + passed_corrected
        stats = {
            "input_reads": total_reads,
            "passed_reads": total_kept,
            "passed_fraction": total_kept / max(1, total_reads),
            "exact_whitelist_matches": passed_exact,
            "exact_match_fraction": passed_exact / max(1, total_reads),
            "corrected_1bp_mismatches": passed_corrected,
            "corrected_fraction": passed_corrected / max(1, total_reads),
            "discarded_invalid_barcodes": discarded_invalid,
            "discarded_short_reads": discarded_short
        }
        
        return stats
