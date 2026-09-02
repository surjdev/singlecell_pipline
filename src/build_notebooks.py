"""
Programmatic Notebook Generator for Smart-seq2 & Downstream Single-Cell Pipeline.
Builds 4 polished, interactive Jupyter notebooks:
1. Raw Data QC & fastp Trimming
2. Genome Alignment (STAR/HISAT2) & Quantification (featureCounts)
3. Single-Cell Count Matrix Assembly & MultiQC
4. Downstream Single-Cell Analysis (Scanpy, QC, PCA, UMAP, Leiden & Markers)
"""

import nbformat as nbf
from pathlib import Path
from src.utils import ensure_dir

def create_notebook_1(nb_dir: Path):
    nb = nbf.v4.new_notebook()
    cells = []
    
    cells.append(nbf.v4.new_markdown_cell(
        "# 🧬 Smart-seq2 Notebook 1: Quality Control & `fastp` Trimming\n\n"
        "This notebook walks through:\n"
        "1. Inspecting raw Smart-seq2 paired-end FASTQ reads\n"
        "2. Automated `FastQC` execution on raw libraries\n"
        "3. High-performance adapter and quality trimming using `fastp` (Nextera adapters, poly-A/G, quality sliding window)\n"
        "4. Before vs. After QC metric comparisons"
    ))
    
    cells.append(nbf.v4.new_code_cell(
        "import os\n"
        "import sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, str(Path('..').resolve()))\n\n"
        "import pandas as pd\n"
        "import matplotlib.pyplot as plt\n"
        "import seaborn as sns\n"
        "from src.utils import load_config\n"
        "from src.sc_qc import SmartSeq2QC\n"
        "from src.sc_preprocess import SmartSeq2Preprocessor\n\n"
        "config = load_config('../config/pipeline_config.yaml')\n"
        "print(f\"Project: {config['project']['name']}\")"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 1. Run Raw FastQC"))
    cells.append(nbf.v4.new_code_cell(
        "qc_engine = SmartSeq2QC(config)\n"
        "qc_engine.run_fastqc_stage('raw')"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 2. Execute `fastp` Paired-End Trimming"))
    cells.append(nbf.v4.new_code_cell(
        "preprocessor = SmartSeq2Preprocessor(config)\n"
        "trim_results = preprocessor.run_all()\n"
        "df_trim = pd.DataFrame(trim_results)\n"
        "df_trim"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 3. Visual Comparison: Read Filtering & Q30 Improvement"))
    cells.append(nbf.v4.new_code_cell(
        "fig, axes = plt.subplots(1, 2, figsize=(12, 4))\n"
        "sns.barplot(data=df_trim, x='cell_id', y='clean_reads', ax=axes[0], color='#2b5c8f')\n"
        "axes[0].set_title('Clean Reads per Cell')\n"
        "axes[0].tick_params(axis='x', rotation=45)\n\n"
        "df_q30 = df_trim.melt(id_vars=['cell_id'], value_vars=['raw_q30_rate', 'clean_q30_rate'], var_name='Stage', value_name='Q30_Rate')\n"
        "sns.barplot(data=df_q30, x='cell_id', y='Q30_Rate', hue='Stage', ax=axes[1], palette='Set2')\n"
        "axes[1].set_title('Q30 Rate Before vs After fastp')\n"
        "axes[1].set_ylim(0, 1.05)\n"
        "axes[1].tick_params(axis='x', rotation=45)\n"
        "plt.tight_layout()\n"
        "plt.show()"
    ))
    
    nb.cells = cells
    out_p = nb_dir / "01_raw_qc_and_trimming.ipynb"
    with open(out_p, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    return out_p

def create_notebook_2(nb_dir: Path):
    nb = nbf.v4.new_notebook()
    cells = []
    
    cells.append(nbf.v4.new_markdown_cell(
        "# 🎯 Smart-seq2 Notebook 2: Splice-Aware Alignment & Quantification\n\n"
        "This notebook walks through:\n"
        "1. Splice-aware reference genome indexing & alignment using **STAR** / **HISAT2**\n"
        "2. Coordinate sorting & BAM indexing via `samtools`\n"
        "3. Gene-level quantification via **Subread `featureCounts`** and **STAR GeneCounts**\n"
        "4. Expression matrix assembly and per-cell alignment rate benchmarking"
    ))
    
    cells.append(nbf.v4.new_code_cell(
        "import os\n"
        "import sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, str(Path('..').resolve()))\n\n"
        "import pandas as pd\n"
        "import matplotlib.pyplot as plt\n"
        "import seaborn as sns\n"
        "from src.utils import load_config\n"
        "from src.aligner import STARAligner, HISAT2Aligner\n"
        "from src.quantifier import FeatureCountsQuantifier, STARCountsParser\n\n"
        "config = load_config('../config/pipeline_config.yaml')"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 1. Run STAR Alignment & BAM Indexing"))
    cells.append(nbf.v4.new_code_cell(
        "star_aligner = STARAligner(config)\n"
        "star_results = star_aligner.run_all()\n"
        "df_star = pd.DataFrame(star_results)\n"
        "df_star"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 2. Run featureCounts Quantification"))
    cells.append(nbf.v4.new_code_cell(
        "fc_quant = FeatureCountsQuantifier(config)\n"
        "count_matrix = fc_quant.run_quantification()\n"
        "count_matrix.head(10)"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 3. Visualizing Uniquely Mapped Reads per Cell"))
    cells.append(nbf.v4.new_code_cell(
        "plt.figure(figsize=(8, 4))\n"
        "sns.barplot(data=df_star, x='cell_id', y='uniquely_mapped_pct', color='#107050')\n"
        "plt.title('STAR Uniquely Mapped Reads Percentage (%)')\n"
        "plt.ylim(0, 100)\n"
        "plt.ylabel('Unique Mapping %')\n"
        "plt.xticks(rotation=45)\n"
        "plt.show()"
    ))
    
    nb.cells = cells
    out_p = nb_dir / "02_alignment_and_quantification.ipynb"
    with open(out_p, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    return out_p

def create_notebook_3(nb_dir: Path):
    nb = nbf.v4.new_notebook()
    cells = []
    
    cells.append(nbf.v4.new_markdown_cell(
        "# 📊 Smart-seq2 Notebook 3: Count Matrix QC & MultiQC\n\n"
        "This notebook walks through:\n"
        "1. Loading the assembled **AnnData** single-cell object (`.h5ad`)\n"
        "2. Calculating single-cell QC metrics (Total Counts, Detected Genes, Mitochondrial Ratio)\n"
        "3. Aggregating all logs into an interactive **MultiQC** dashboard"
    ))
    
    cells.append(nbf.v4.new_code_cell(
        "import os\n"
        "import sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, str(Path('..').resolve()))\n\n"
        "import anndata as ad\n"
        "import pandas as pd\n"
        "import matplotlib.pyplot as plt\n"
        "import seaborn as sns\n"
        "from src.utils import load_config\n"
        "from src.sc_qc import SmartSeq2QC\n\n"
        "config = load_config('../config/pipeline_config.yaml')\n"
        "adata = ad.read_h5ad('../data/counts/gene_cell_count_matrix.h5ad')\n"
        "adata"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 1. Single-Cell Per-Cell QC Summary"))
    cells.append(nbf.v4.new_code_cell(
        "adata.obs"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 2. QC Distributions (Total Counts & Detected Genes)"))
    cells.append(nbf.v4.new_code_cell(
        "fig, axes = plt.subplots(1, 2, figsize=(10, 4))\n"
        "sns.scatterplot(data=adata.obs, x='total_counts', y='n_genes_by_counts', ax=axes[0], color='#9b2226', s=70)\n"
        "axes[0].set_title('Total Counts vs. Detected Genes')\n"
        "axes[0].set_xlabel('Total Read Counts')\n"
        "axes[0].set_ylabel('Number of Detected Genes')\n\n"
        "sns.boxplot(data=adata.obs[['total_counts', 'n_genes_by_counts']], ax=axes[1], palette='Pastel1')\n"
        "axes[1].set_title('Distribution of Cell Metrics')\n"
        "plt.tight_layout()\n"
        "plt.show()"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 3. Generate MultiQC Dashboard"))
    cells.append(nbf.v4.new_code_cell(
        "qc_engine = SmartSeq2QC(config)\n"
        "qc_engine.run_multiqc()\n"
        "print(\"MultiQC report compiled in reports/multiqc/multiqc_report.html\")"
    ))
    
    nb.cells = cells
    out_p = nb_dir / "03_single_cell_matrix_and_multiqc.ipynb"
    with open(out_p, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    return out_p

def create_notebook_4(nb_dir: Path):
    nb = nbf.v4.new_notebook()
    cells = []
    
    cells.append(nbf.v4.new_markdown_cell(
        "# 🔬 Notebook 4: Downstream Single-Cell Analysis (Scanpy, UMAP, Leiden & Markers)\n\n"
        "This notebook walks through the complete downstream single-cell analysis workflow:\n"
        "1. **Quality Control & Filtering**: Cell & gene filtering, mitochondrial percentage thresholds\n"
        "2. **Normalization & Log1p**: Library size scaling & variance stabilization\n"
        "3. **Highly Variable Genes (HVGs)**: Feature selection\n"
        "4. **Dimensionality Reduction**: Principal Component Analysis (PCA) & Neighborhood Graph\n"
        "5. **Non-Linear Embedding**: UMAP visualization\n"
        "6. **Cell Clustering**: Leiden community detection algorithm\n"
        "7. **Marker Gene Discovery**: Cluster-specific differential expression (Wilcoxon rank-sum test)\n"
        "8. **Visualization**: Publication-quality UMAPs, DotPlots, and Heatmaps"
    ))
    
    cells.append(nbf.v4.new_code_cell(
        "import os\n"
        "import sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, str(Path('..').resolve()))\n\n"
        "import scanpy as sc\n"
        "import anndata as ad\n"
        "import pandas as pd\n"
        "import matplotlib.pyplot as plt\n"
        "import seaborn as sns\n"
        "from src.utils import load_config\n"
        "from src.downstream import SingleCellDownstream\n\n"
        "sc.set_figure_params(dpi=120, facecolor='white', frameon=True)\n"
        "config = load_config('../config/pipeline_config.yaml')"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 1. Load Count Matrix & Compute QC Metrics"))
    cells.append(nbf.v4.new_code_cell(
        "ds = SingleCellDownstream(config)\n"
        "adata = ds.load_anndata()\n"
        "adata = ds.run_qc_and_filtering(adata)\n"
        "adata"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 2. Normalization & HVG Selection"))
    cells.append(nbf.v4.new_code_cell(
        "adata = ds.run_normalization_and_hvg(adata)\n"
        "print(f\"Highly variable genes: {adata.var['highly_variable'].sum()} / {adata.n_vars}\")"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 3. PCA, Neighborhood Graph, UMAP & Leiden Clustering"))
    cells.append(nbf.v4.new_code_cell(
        "adata = ds.run_clustering_and_embeddings(adata)\n"
        "sc.pl.umap(adata, color=['leiden', 'total_counts', 'n_genes_by_counts'], ncols=3)"
    ))
    
    cells.append(nbf.v4.new_markdown_cell("### 4. Cluster Marker Gene Discovery (Differential Expression)"))
    cells.append(nbf.v4.new_code_cell(
        "df_markers = ds.run_marker_discovery(adata)\n"
        "if not df_markers.empty:\n"
        "    display(df_markers.head(10))\n"
        "    top_genes = list(df_markers.groupby('cluster')['gene'].first().unique())\n"
        "    sc.pl.dotplot(adata, var_names=top_genes, groupby='leiden')"
    ))
    
    nb.cells = cells
    out_p = nb_dir / "04_downstream_single_cell_analysis.ipynb"
    with open(out_p, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    return out_p

def main():
    nb_dir = ensure_dir("notebooks")
    create_notebook_1(nb_dir)
    create_notebook_2(nb_dir)
    create_notebook_3(nb_dir)
    create_notebook_4(nb_dir)
    print("✔ Successfully built all 4 Smart-seq2 & Downstream notebooks!")

if __name__ == "__main__":
    main()
