"""
Single-Cell FASTQ Preprocessing and Quality Control Pipeline CLI.
Orchestrates raw QC, quality trimming, barcode error correction, clean QC,
and MultiQC aggregated reporting.
"""

import sys
import json
import time
import argparse
from pathlib import Path
from rich.console import Console
from rich.panel import Panel

from src.utils import load_config, ensure_dir, run_cmd, print_summary_table
from src.sc_qc import SingleCellQC
from src.sc_preprocess import SingleCellPreprocessor

console = Console()

def run_qc_stage(config: dict, stage: str = "raw") -> dict:
    """Run FastQC and single-cell QC analysis for a given stage (raw or clean)."""
    paths = config["paths"]
    chem = config.get("chemistry", {})
    cb_len = chem.get("r1_structure", {}).get("cell_barcode_len", 16)
    umi_len = chem.get("r1_structure", {}).get("umi_len", 12)
    sample_id = config.get("project", {}).get("sample_id", "sample")
    
    if stage == "raw":
        r1 = paths["sample_r1"]
        r2 = paths["sample_r2"]
        out_dir = paths["qc_raw_dir"]
    else:
        r1 = paths["clean_r1"]
        r2 = paths["clean_r2"]
        out_dir = paths["qc_clean_dir"]
        
    console.print(Panel.fit(f"[bold cyan]Running QC for Stage: [magenta]{stage.upper()}[/magenta][/bold cyan]"))
    
    if not Path(r1).exists() or not Path(r2).exists():
        raise FileNotFoundError(f"Input FASTQ files for stage '{stage}' not found: {r1}, {r2}")
        
    qc = SingleCellQC(
        r1_path=r1,
        r2_path=r2,
        cb_len=cb_len,
        umi_len=umi_len,
        whitelist_path=paths.get("whitelist_file"),
        sample_id=f"{sample_id}_{stage}"
    )
    
    # 1. Run FastQC
    qc.run_fastqc(out_dir=out_dir)
    
    # 2. Run Single Cell metric computation
    results = qc.analyze_fastq()
    
    # 3. Save report JSON
    json_path = Path(out_dir) / f"{sample_id}_{stage}_sc_qc.json"
    qc.save_report(results, json_path)
    
    # Print summary
    print_summary_table(f"Single-Cell QC Metrics ({stage.upper()})", results["summary_metrics"])
    return results

def run_preprocessing_stage(config: dict) -> dict:
    """Run fastp trimming and barcode error correction."""
    paths = config["paths"]
    sample_id = config.get("project", {}).get("sample_id", "sample")
    
    raw_r1 = Path(paths["sample_r1"])
    raw_r2 = Path(paths["sample_r2"])
    
    if not raw_r1.exists() or not raw_r2.exists():
        raise FileNotFoundError(f"Raw FASTQ files not found: {raw_r1}, {raw_r2}")
        
    console.print(Panel.fit("[bold green]Starting Preprocessing & Quality Trimming[/bold green]"))
    
    preprocessor = SingleCellPreprocessor(config)
    
    # Temp trimmed paths
    temp_dir = ensure_dir(Path("data/temp"))
    trimmed_r1 = temp_dir / f"{sample_id}_trimmed_R1.fastq.gz"
    trimmed_r2 = temp_dir / f"{sample_id}_trimmed_R2.fastq.gz"
    fp_html = Path(paths["reports_dir"]) / "fastp" / f"{sample_id}_fastp.html"
    fp_json = Path(paths["reports_dir"]) / "fastp" / f"{sample_id}_fastp.json"
    
    # 1. Fastp trimming
    fp_stats = preprocessor.run_fastp_trimming(
        r1_in=raw_r1,
        r2_in=raw_r2,
        r1_out=trimmed_r1,
        r2_out=trimmed_r2,
        report_html=fp_html,
        report_json=fp_json
    )
    
    # 2. Barcode error-correction and final clean FASTQ output
    clean_r1 = Path(paths["clean_r1"])
    clean_r2 = Path(paths["clean_r2"])
    extracted_out = Path(paths.get("extracted_fastq")) if paths.get("extracted_fastq") else None
    
    bc_stats = preprocessor.filter_and_correct_barcodes(
        r1_in=trimmed_r1,
        r2_in=trimmed_r2,
        r1_out=clean_r1,
        r2_out=clean_r2,
        extracted_out=extracted_out,
        whitelist_file=paths.get("whitelist_file")
    )
    
    # Save combined preprocessing stats JSON
    prep_stats_file = Path(paths["reports_dir"]) / f"{sample_id}_preprocessing_summary.json"
    prep_stats = {
        "fastp_summary": {
            "before_filtering_reads": fp_stats.get("summary", {}).get("before_filtering", {}).get("total_reads", 0),
            "after_filtering_reads": fp_stats.get("summary", {}).get("after_filtering", {}).get("total_reads", 0),
            "filtering_passed_reads_fraction": fp_stats.get("filtering_result", {}).get("passed_filter_reads", 0) / max(1, fp_stats.get("summary", {}).get("before_filtering", {}).get("total_reads", 1)),
            "q30_rate_before": fp_stats.get("summary", {}).get("before_filtering", {}).get("q30_rate", 0),
            "q30_rate_after": fp_stats.get("summary", {}).get("after_filtering", {}).get("q30_rate", 0),
            "adapter_trimmed_reads": fp_stats.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
            "polyx_trimmed_reads": fp_stats.get("polyx_trimming", {}).get("total_polyx_trimmed_reads", 0)
        },
        "barcode_filtering_summary": bc_stats
    }
    
    with open(prep_stats_file, "w") as f:
        json.dump(prep_stats, f, indent=2)
        
    print_summary_table("Barcode Filtering & Error Correction", bc_stats)
    
    # Clean up temp
    if trimmed_r1.exists():
        trimmed_r1.unlink()
    if trimmed_r2.exists():
        trimmed_r2.unlink()
    if temp_dir.exists() and not list(temp_dir.iterdir()):
        temp_dir.rmdir()
        
    return prep_stats

def run_multiqc_stage(config: dict) -> Path:
    """Run MultiQC to aggregate FastQC and Fastp reports into a single interactive HTML report."""
    paths = config["paths"]
    reports_dir = Path(paths["reports_dir"])
    multiqc_out = ensure_dir(paths["multiqc_dir"])
    
    console.print(Panel.fit("[bold magenta]Generating Aggregated MultiQC Dashboard[/bold magenta]"))
    
    cmd = [
        "multiqc",
        str(reports_dir),
        "-o", str(multiqc_out),
        "--title", f"Single-Cell QC Dashboard - {config.get('project', {}).get('name', 'scRNA')}",
        "--filename", "single_cell_multiqc_report.html",
        "--force"
    ]
    
    run_cmd(cmd, desc="Running MultiQC")
    report_file = multiqc_out / "single_cell_multiqc_report.html"
    console.print(f"[bold green]✔ MultiQC Report Generated:[/bold green] [cyan]{report_file}[/cyan]")
    return report_file

def run_full_pipeline(config_path: str = "config/pipeline_config.yaml"):
    """Execute end-to-end Single-Cell RAW to Clean pipeline."""
    start_time = time.time()
    config = load_config(config_path)
    
    console.print(Panel.fit("[bold blue]Starting Full Single-Cell RAW-to-Clean FASTQ Pipeline[/bold blue]"))
    
    # Stage 1: Raw QC
    raw_qc = run_qc_stage(config, stage="raw")
    
    # Stage 2: Preprocessing & Trimming
    prep_stats = run_preprocessing_stage(config)
    
    # Stage 3: Clean QC
    clean_qc = run_qc_stage(config, stage="clean")
    
    # Stage 4: MultiQC
    multiqc_report = run_multiqc_stage(config)
    
    elapsed = time.time() - start_time
    
    # Final Comparison Summary
    raw_m = raw_qc["summary_metrics"]
    clean_m = clean_qc["summary_metrics"]
    
    comparison = {
        "Raw Reads": raw_m["total_reads"],
        "Clean Reads": clean_m["total_reads"],
        "Clean Yield Rate": clean_m["total_reads"] / max(1, raw_m["total_reads"]),
        "Raw R2 Q30 Fraction": raw_m["r2_q30_fraction"],
        "Clean R2 Q30 Fraction": clean_m["r2_q30_fraction"],
        "Raw Valid Barcodes": raw_m.get("valid_barcode_fraction", 0.0),
        "Clean Valid Barcodes": clean_m.get("valid_barcode_fraction", 0.0),
        "Execution Time (s)": f"{elapsed:.2f}s"
    }
    
    print_summary_table("End-to-End Pipeline Summary (Before vs After)", comparison)
    console.print("[bold green]✔ Pipeline finished successfully![/bold green]")

def main():
    parser = argparse.ArgumentParser(description="Single-cell FASTQ Preprocessing & QC Pipeline")
    parser.add_argument("--config", default="config/pipeline_config.yaml", help="Path to pipeline configuration YAML")
    subparsers = parser.add_subparsers(dest="command", help="Pipeline subcommands")
    
    # Subcommand: qc
    qc_parser = subparsers.add_parser("qc", help="Run QC on raw or clean FASTQs")
    qc_parser.add_argument("--stage", choices=["raw", "clean"], default="raw", help="QC stage to run")
    
    # Subcommand: preprocess
    subparsers.add_parser("preprocess", help="Run fastp trimming and barcode error correction")
    
    # Subcommand: multiqc
    subparsers.add_parser("multiqc", help="Run MultiQC report aggregation")
    
    # Subcommand: full
    subparsers.add_parser("full", help="Run complete end-to-end pipeline")
    
    args = parser.parse_args()
    
    if not args.command or args.command == "full":
        run_full_pipeline(args.config)
    elif args.command == "qc":
        cfg = load_config(args.config)
        run_qc_stage(cfg, stage=args.stage)
    elif args.command == "preprocess":
        cfg = load_config(args.config)
        run_preprocessing_stage(cfg)
    elif args.command == "multiqc":
        cfg = load_config(args.config)
        run_multiqc_stage(cfg)

if __name__ == "__main__":
    main()
