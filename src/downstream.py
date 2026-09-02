"""
Downstream Single-Cell Analysis Engine.
Implements complete single-cell workflow using Scanpy:
- Cell & Gene QC filtering
- Depth normalization & Log1p transformation
- Highly Variable Gene (HVG) selection
- PCA dimensionality reduction & Neighbor Graph construction
- Non-linear UMAP embedding & Leiden graph clustering
- Marker Gene Discovery (Differential Expression per cluster)
- Publication-quality plot generation and AnnData export.
"""

import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Any, Optional
from rich.console import Console
from rich.table import Table

from src.utils import load_config, ensure_dir

console = Console()

class SingleCellDownstream:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.paths = config["paths"]
        self.ds_cfg = config.get("downstream", {})
        
        self.input_h5ad = Path(self.paths.get("input_h5ad", "data/counts/gene_cell_count_matrix.h5ad"))
        self.out_h5ad = Path(self.paths.get("processed_h5ad", "data/counts/processed_adata.h5ad"))
        self.out_dir = ensure_dir(self.paths.get("downstream_dir", "reports/downstream"))
        
        # Configure Scanpy plot settings
        sc.settings.verbosity = 1
        sc.settings.figdir = str(self.out_dir)
        sc.set_figure_params(dpi=150, facecolor="white", frameon=True)

    def load_anndata(self) -> ad.AnnData:
        """Load the raw count matrix AnnData object."""
        if not self.input_h5ad.exists():
            raise FileNotFoundError(f"Count matrix AnnData not found at {self.input_h5ad}. Run quantification step first.")
            
        adata = sc.read_h5ad(self.input_h5ad)
        console.print(f"[bold cyan]Loaded count matrix:[/bold cyan] {adata.n_obs} cells × {adata.n_vars} genes")
        return adata

    def run_qc_and_filtering(self, adata: ad.AnnData) -> ad.AnnData:
        """Calculate QC metrics and apply cell/gene filters."""
        filt_cfg = self.ds_cfg.get("filtering", {})
        
        # Calculate QC metrics
        is_mito = adata.var_names.str.contains("MT-", case=False) | adata.var_names.str.contains("MT_", case=False)
        adata.var["mt"] = is_mito
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        
        # Save QC violin plot
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
        sns.violinplot(y=adata.obs["total_counts"], ax=axes[0], color="#2b5c8f")
        axes[0].set_title("Total Counts per Cell")
        sns.violinplot(y=adata.obs["n_genes_by_counts"], ax=axes[1], color="#107050")
        axes[1].set_title("Detected Genes per Cell")
        sns.violinplot(y=adata.obs["pct_counts_mt"], ax=axes[2], color="#9b2226")
        axes[2].set_title("Mitochondrial %")
        plt.tight_layout()
        qc_plot_path = self.out_dir / "qc_violin.png"
        fig.savefig(qc_plot_path)
        plt.close(fig)
        
        # Filtering thresholds
        min_genes = filt_cfg.get("min_genes_per_cell", 3)
        min_cells = filt_cfg.get("min_cells_per_gene", 1)
        max_mito = filt_cfg.get("max_pct_mito", 25.0)
        
        n_cells_before = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata = adata[adata.obs["pct_counts_mt"] <= max_mito].copy()
        
        console.print(f"[green]✔ QC filtering complete:[/green] Retained {adata.n_obs}/{n_cells_before} cells and {adata.n_vars} genes")
        return adata

    def run_normalization_and_hvg(self, adata: ad.AnnData) -> ad.AnnData:
        """Perform library size normalization, log1p transformation, and HVG selection."""
        norm_cfg = self.ds_cfg.get("normalization", {})
        hvg_cfg = self.ds_cfg.get("feature_selection", {})
        
        # Store raw counts in layer
        adata.layers["raw_counts"] = adata.X.copy()
        
        # Normalize total counts
        target_sum = norm_cfg.get("target_sum", 1e4)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        
        if norm_cfg.get("log1p", True):
            sc.pp.log1p(adata)
            
        # Store normalized data in .raw
        adata.raw = adata
        
        # Highly Variable Genes
        n_top = min(hvg_cfg.get("n_top_genes", 2000), adata.n_vars)
        if adata.n_vars > 10:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=hvg_cfg.get("flavor", "seurat"))
        else:
            adata.var["highly_variable"] = True
            
        console.print(f"[green]✔ Normalization & HVG selection complete[/green]")
        return adata

    def run_clustering_and_embeddings(self, adata: ad.AnnData) -> ad.AnnData:
        """Perform PCA, neighborhood graph construction, UMAP, and Leiden clustering."""
        dim_cfg = self.ds_cfg.get("dim_reduction", {})
        clust_cfg = self.ds_cfg.get("clustering", {})
        
        # Scale data
        sc.pp.scale(adata, max_value=10)
        
        # PCA
        n_pcs = min(dim_cfg.get("n_pcs", 10), adata.n_obs - 1, adata.n_vars - 1)
        n_pcs = max(2, n_pcs)
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
        
        # Neighbors
        n_neighbors = min(dim_cfg.get("n_neighbors", 5), adata.n_obs - 1)
        n_neighbors = max(2, n_neighbors)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        
        # UMAP
        sc.tl.umap(adata)
        
        # Leiden clustering
        resolution = clust_cfg.get("resolution", 0.5)
        sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
        
        # Save UMAP plot
        umap_plot_path = self.out_dir / "umap_clusters.png"
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color="leiden", title="Single-Cell Clusters (Leiden)", ax=ax, show=False)
        fig.savefig(umap_plot_path, bbox_inches="tight")
        plt.close(fig)
        
        n_clusters = len(adata.obs["leiden"].unique())
        console.print(f"[green]✔ Clustering complete:[/green] Identified [bold magenta]{n_clusters}[/bold magenta] Leiden cluster(s)")
        return adata

    def run_marker_discovery(self, adata: ad.AnnData) -> pd.DataFrame:
        """Identify cluster marker genes using differential expression."""
        marker_cfg = self.ds_cfg.get("marker_discovery", {})
        method = marker_cfg.get("method", "wilcoxon")
        n_genes = marker_cfg.get("n_genes", 10)
        
        if len(adata.obs["leiden"].unique()) > 1:
            sc.tl.rank_genes_groups(adata, groupby="leiden", method=method, use_raw=True)
            
            # Export marker table
            result = adata.uns["rank_genes_groups"]
            groups = result["names"].dtype.names
            
            records = []
            for group in groups:
                for rank in range(min(n_genes, len(result["names"][group]))):
                    records.append({
                        "cluster": group,
                        "rank": rank + 1,
                        "gene": result["names"][group][rank],
                        "score": result["scores"][group][rank],
                        "logfoldchange": result["logfoldchanges"][group][rank],
                        "pvals_adj": result["pvals_adj"][group][rank]
                    })
                    
            df_markers = pd.DataFrame(records)
            markers_csv = self.out_dir / "cluster_markers.csv"
            df_markers.to_csv(markers_csv, index=False)
            
            # Save DotPlot if multiple clusters
            try:
                top_genes = list(df_markers.groupby("cluster")["gene"].first().unique())
                fig_dp = sc.pl.dotplot(adata, var_names=top_genes, groupby="leiden", return_fig=True, show=False)
                dp_path = self.out_dir / "marker_dotplot.png"
                fig_dp.savefig(dp_path, bbox_inches="tight")
            except Exception:
                pass
                
            console.print(f"[green]✔ Top marker genes discovered and saved to:[/green] {markers_csv}")
            return df_markers
        else:
            console.print("[yellow]Single cluster detected; skipping differential marker discovery.[/yellow]")
            return pd.DataFrame()

    def run_all(self) -> ad.AnnData:
        """Execute the full downstream single-cell pipeline."""
        console.print(f"[bold cyan]Starting Downstream Single-Cell Analysis...[/bold cyan]")
        adata = self.load_anndata()
        adata = self.run_qc_and_filtering(adata)
        adata = self.run_normalization_and_hvg(adata)
        adata = self.run_clustering_and_embeddings(adata)
        df_markers = self.run_marker_discovery(adata)
        
        # Save processed AnnData
        adata.write_h5ad(self.out_h5ad)
        console.print(f"[bold green]✔ Processed AnnData saved:[/bold green] {self.out_h5ad}")
        
        # Display Summary Table
        self._print_summary(adata, df_markers)
        return adata

    def _print_summary(self, adata: ad.AnnData, df_markers: pd.DataFrame):
        table = Table(title="Downstream Single-Cell Analysis Summary", header_style="bold magenta", border_style="cyan")
        table.add_column("Metric", style="bold white")
        table.add_column("Value", style="green", justify="right")
        
        table.add_row("Total Cells Analyzed", f"{adata.n_obs:,}")
        table.add_row("Total Genes Expressed", f"{adata.n_vars:,}")
        table.add_row("Leiden Clusters Identified", f"{len(adata.obs['leiden'].unique()):,}")
        table.add_row("Top Marker Genes Exported", f"{len(df_markers):,}")
        table.add_row("UMAP Plot Saved", str(self.out_dir / "umap_clusters.png"))
        table.add_row("QC Violin Plot Saved", str(self.out_dir / "qc_violin.png"))
        console.print(table)

def main():
    config = load_config()
    ds = SingleCellDownstream(config)
    ds.run_all()

if __name__ == "__main__":
    main()
