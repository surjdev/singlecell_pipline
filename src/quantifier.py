"""
Transcriptomic & Single-Cell Quantification Engine.
Wraps featureCounts (Subread) and STAR ReadsPerGene count parsers.
Aggregates per-cell count profiles into a unified Gene x Cell Expression Matrix
(CSV, TSV, and AnnData .h5ad) and computes per-cell QC statistics.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional
import anndata as ad
from rich.console import Console
from rich.table import Table

from src.utils import run_cmd, ensure_dir, load_config

console = Console()

class FeatureCountsQuantifier:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.quant_cfg = config.get("quantification", {})
        self.fc_cfg = self.quant_cfg.get("featurecounts", {})
        self.threads = self.quant_cfg.get("threads", 4)
        
        self.aligned_dir = Path(self.paths.get("aligned_dir", "data/aligned"))
        self.counts_dir = ensure_dir(self.paths.get("counts_dir", "data/counts"))
        self.gtf_file = Path(self.paths["gtf_file"])

    def run_quantification(self, aligner_pref: Optional[str] = None) -> pd.DataFrame:
        """Run featureCounts on aligned BAM files and assemble the count matrix."""
        pref = aligner_pref or self.config.get("alignment", {}).get("default_aligner", "star").lower()
        
        if pref == "star":
            bam_files = sorted(list(self.aligned_dir.glob("*_Aligned.sortedByCoord.out.bam")))
        elif pref == "hisat2":
            bam_files = sorted(list(self.aligned_dir.glob("*_hisat2_sorted.bam")))
        else:
            bam_files = []
            
        if not bam_files:
            bam_files = sorted(list(self.aligned_dir.glob("*.bam")))
            bam_files = [b for b in bam_files if not b.name.endswith(".tmp.bam")]
        
        if not bam_files:
            raise FileNotFoundError(f"No BAM files found in {self.aligned_dir} for featureCounts.")
            
        out_raw_counts = self.counts_dir / "featureCounts_raw.txt"
        summary_file = self.counts_dir / "featureCounts_raw.txt.summary"
        
        cmd = [
            "featureCounts",
            "-a", str(self.gtf_file),
            "-o", str(out_raw_counts),
            "-T", str(self.threads),
            "-t", str(self.fc_cfg.get("feature_type", "exon")),
            "-g", str(self.fc_cfg.get("attribute_type", "gene_id")),
            "-s", str(self.fc_cfg.get("strand_specificity", 0))
        ]
        
        if self.fc_cfg.get("is_paired_end", True):
            cmd.append("-p")
        if self.fc_cfg.get("count_read_pairs", True):
            cmd.append("--countReadPairs")
            
        extra_args = self.fc_cfg.get("extra_args")
        if extra_args:
            cmd.extend(extra_args.split())
            
        cmd.extend([str(b) for b in bam_files])
        run_cmd(cmd, desc="Running featureCounts across cell BAMs")
        
        # Parse featureCounts table into clean DataFrame
        df_raw = pd.read_csv(out_raw_counts, sep="\t", comment="#")
        gene_ids = df_raw["Geneid"]
        count_cols = df_raw.columns[6:]
        
        matrix_df = df_raw[count_cols].copy()
        matrix_df.index = gene_ids
        
        # Clean column names from full paths to cell IDs
        clean_col_names = []
        for col in count_cols:
            col_path = Path(col)
            c_name = col_path.name.replace("_Aligned.sortedByCoord.out.bam", "").replace("_hisat2_sorted.bam", "").replace(".bam", "")
            clean_col_names.append(c_name)
        matrix_df.columns = clean_col_names
        
        self._export_matrix(matrix_df, "featurecounts")
        return matrix_df

    def _export_matrix(self, df: pd.DataFrame, source: str):
        csv_p = self.counts_dir / "gene_cell_count_matrix.csv"
        tsv_p = self.counts_dir / "gene_cell_count_matrix.tsv"
        h5ad_p = self.counts_dir / "gene_cell_count_matrix.h5ad"
        
        df.to_csv(csv_p)
        df.to_csv(tsv_p, sep="\t")
        
        # Build AnnData object (Cells x Genes)
        adata = ad.AnnData(X=df.T.values, obs=pd.DataFrame(index=df.columns.astype(str)), var=pd.DataFrame(index=df.index.astype(str)))
        adata.obs_names_make_unique()
        
        # Compute cell-level QC
        adata.obs["total_counts"] = np.asarray(np.sum(adata.X, axis=1)).flatten()
        adata.obs["n_genes_by_counts"] = np.asarray(np.sum(adata.X > 0, axis=1)).flatten()
        
        # Check mitochondrial genes
        is_mito = adata.var_names.str.contains("MT-", case=False) | adata.var_names.str.contains("MT_", case=False)
        if np.any(is_mito):
            mito_counts = np.asarray(np.sum(adata.X[:, is_mito], axis=1)).flatten()
            adata.obs["pct_counts_mt"] = (mito_counts / np.maximum(1, adata.obs["total_counts"].values)) * 100.0
        else:
            adata.obs["pct_counts_mt"] = 0.0
            
        adata.write_h5ad(h5ad_p)
        
        console.print(f"[bold green]✔ Saved Count Matrix ({df.shape[0]} genes × {df.shape[1]} cells):[/bold green]")
        console.print(f"  - CSV: {csv_p}")
        console.print(f"  - TSV: {tsv_p}")
        console.print(f"  - AnnData: {h5ad_p}")
        
        # Display Single-Cell QC table
        self._print_sc_qc(adata)

    def _print_sc_qc(self, adata: ad.AnnData):
        table = Table(title="Single-Cell Expression QC Metrics", header_style="bold magenta", border_style="cyan")
        table.add_column("Cell ID", style="bold white")
        table.add_column("Total Counts", justify="right", style="green")
        table.add_column("Detected Genes", justify="right", style="cyan")
        table.add_column("Mito %", justify="right", style="yellow")
        
        for i, cell_id in enumerate(adata.obs_names):
            tot = float(adata.obs["total_counts"].iloc[i])
            ngen = int(adata.obs["n_genes_by_counts"].iloc[i])
            pct_mt = float(adata.obs["pct_counts_mt"].iloc[i])
            table.add_row(cell_id, f"{int(tot):,}", f"{ngen:,}", f"{pct_mt:.2f}%")
        console.print(table)



class STARCountsParser:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.aligned_dir = Path(self.paths.get("aligned_dir", "data/aligned"))
        self.counts_dir = ensure_dir(self.paths.get("counts_dir", "data/counts"))
        
    def run_quantification(self) -> pd.DataFrame:
        """Parse ReadsPerGene.out.tab files generated by STAR."""
        tab_files = sorted(list(self.aligned_dir.glob("*_ReadsPerGene.out.tab")))
        if not tab_files:
            raise FileNotFoundError(f"No *_ReadsPerGene.out.tab files found in {self.aligned_dir}. Ensure STAR was run with --quantMode GeneCounts.")
            
        cell_dict = {}
        for tab in tab_files:
            cell_id = tab.name.replace("_ReadsPerGene.out.tab", "")
            # ReadsPerGene format:
            # line 0-3: N_unmapped, N_multimapping, N_noFeature, N_ambiguous
            # col 0: gene_id, col 1: unstranded, col 2: 1st read strand, col 3: 2nd read strand
            df_cell = pd.read_csv(tab, sep="\t", header=None, skiprows=4, index_col=0)
            cell_dict[cell_id] = df_cell[1] # Unstranded counts
            
        matrix_df = pd.DataFrame(cell_dict)
        matrix_df.index.name = "Geneid"
        
        fc = FeatureCountsQuantifier(self.config)
        fc._export_matrix(matrix_df, "star_counts")
        return matrix_df

def main():
    config = load_config()
    method = config.get("quantification", {}).get("default_method", "featurecounts").lower()
    if method == "star":
        quant = STARCountsParser(config)
    else:
        quant = FeatureCountsQuantifier(config)
    quant.run_quantification()

if __name__ == "__main__":
    main()
