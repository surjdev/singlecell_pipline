"""
Smart-seq2 & Full-Length scRNA-seq Pipeline Orchestrator.
Unified CLI workflow supporting FastQC, fastp, STAR, HISAT2,
featureCounts, and interactive MultiQC reporting.
"""

import sys
import click
from pathlib import Path
from rich.console import Console
from rich.panel import Panel

from src.utils import load_config, ensure_dir
from src.sc_preprocess import SmartSeq2Preprocessor
from src.aligner import STARAligner, HISAT2Aligner
from src.quantifier import FeatureCountsQuantifier, STARCountsParser
from src.sc_qc import SmartSeq2QC
from src.downstream import SingleCellDownstream

console = Console()

@click.group()
@click.option("--config", default="config/pipeline_config.yaml", help="Path to YAML configuration file.")
@click.pass_context
def cli(ctx, config):
    """Smart-seq2 & Transcriptomic Single-Cell Processing CLI."""
    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config
    ctx.obj["config"] = load_config(config)

@cli.command()
@click.option("--stage", type=click.Choice(["raw", "clean"]), default="raw", help="FASTQ stage to analyze.")
@click.pass_context
def qc(ctx, stage):
    """Run FastQC quality control analysis."""
    config = ctx.obj["config"]
    qc_engine = SmartSeq2QC(config)
    qc_engine.run_fastqc_stage(stage)

@cli.command()
@click.pass_context
def trim(ctx):
    """Run fastp paired-end adapter and quality trimming."""
    config = ctx.obj["config"]
    preprocessor = SmartSeq2Preprocessor(config)
    preprocessor.run_all()

@cli.command()
@click.option("--aligner", type=click.Choice(["star", "hisat2"]), default=None, help="Aligner to use (STAR or HISAT2).")
@click.pass_context
def align(ctx, aligner):
    """Run splice-aware genome alignment with STAR or HISAT2."""
    config = ctx.obj["config"]
    chosen = aligner or config.get("alignment", {}).get("default_aligner", "star").lower()
    if chosen == "hisat2":
        aln = HISAT2Aligner(config)
    else:
        aln = STARAligner(config)
    aln.run_all()

@cli.command()
@click.option("--method", type=click.Choice(["featurecounts", "star"]), default=None, help="Quantification method.")
@click.pass_context
def quant(ctx, method):
    """Run transcriptomic / single-cell gene quantification."""
    config = ctx.obj["config"]
    chosen = method or config.get("quantification", {}).get("default_method", "featurecounts").lower()
    if chosen == "star":
        q = STARCountsParser(config)
    else:
        q = FeatureCountsQuantifier(config)
    q.run_quantification()

@cli.command()
@click.pass_context
def multiqc(ctx):
    """Aggregate all QC, trimming, alignment, and quantification logs into MultiQC."""
    config = ctx.obj["config"]
    qc_engine = SmartSeq2QC(config)
    qc_engine.run_multiqc()

@cli.command()
@click.pass_context
def downstream(ctx):
    """Run downstream single-cell analysis (Scanpy, QC, Normalization, PCA, UMAP, Leiden, Markers)."""
    config = ctx.obj["config"]
    ds = SingleCellDownstream(config)
    ds.run_all()

@cli.command()
@click.option("--aligner", type=click.Choice(["star", "hisat2"]), default="star", help="Aligner to use.")
@click.option("--quant-method", type=click.Choice(["featurecounts", "star"]), default="featurecounts", help="Quantification method.")
@click.option("--run-downstream/--no-downstream", default=True, help="Whether to run downstream Scanpy analysis.")
@click.pass_context
def full(ctx, aligner, quant_method, run_downstream):
    """Execute the complete end-to-end Smart-seq2 pipeline."""
    config = ctx.obj["config"]
    console.print(Panel.fit(
        "[bold cyan]Smart-seq2 & Full-Length scRNA-seq Pipeline Execution[/bold cyan]\n"
        f"Project: [yellow]{config.get('project', {}).get('name')}[/yellow] | Aligner: [green]{aligner}[/green] | Quant: [green]{quant_method}[/green]",
        border_style="magenta"
    ))
    
    # 1. Raw FastQC
    qc_engine = SmartSeq2QC(config)
    qc_engine.run_fastqc_stage("raw")
    
    # 2. fastp Trimming
    preprocessor = SmartSeq2Preprocessor(config)
    preprocessor.run_all()
    
    # 3. Clean FastQC
    qc_engine.run_fastqc_stage("clean")
    
    # 4. Alignment
    if aligner == "hisat2":
        aln = HISAT2Aligner(config)
    else:
        aln = STARAligner(config)
    aln.run_all()
    
    # 5. Quantification
    if quant_method == "star":
        q = STARCountsParser(config)
    else:
        q = FeatureCountsQuantifier(config)
    q.run_quantification()
    
    # 6. MultiQC
    qc_engine.run_multiqc()
    
    # 7. Downstream Single-Cell Analysis
    if run_downstream:
        ds = SingleCellDownstream(config)
        ds.run_all()
        
    console.print("[bold green]✨ Full Smart-seq2 pipeline completed successfully![/bold green]")


def main():
    cli(obj={})

if __name__ == "__main__":
    main()
