# 🧬 Smart-seq2 Single-Cell & Transcriptomic Pipeline (Upstream & Downstream)

A production-grade, reproducible pipeline and interactive Jupyter analysis environment for processing **Smart-seq2 / Plate-based Single-Cell RNA-seq and Transcriptomic data** from raw paired-end FASTQ reads to alignments, gene quantification count matrices, **MultiQC** reporting, and complete **Downstream Single-Cell Analysis (Scanpy, UMAP, Leiden Clustering, Marker Gene Discovery, and GSEA)**.

Managed via **[Pixi](https://pixi.sh/)** with native `bioconda` and `conda-forge` environments.

---

## 🚀 Key Features

### 1. Upstream Processing & Alignment
- **FastQC & MultiQC**: Automated raw and clean read quality assessments, per-tile Phred scores, adapter content, and alignment stats.
- **fastp Trimming**: High-performance Nextera transposase adapter (`CTGTCTCTTATACACATCT`) removal, ISPCR oligo trimming, poly-A/poly-G tail removal, and sliding-window quality filtering.
- **STAR & HISAT2 Alignment**: Splice-aware reference genome indexing, paired-end mapping, coordinate sorting, and BAM indexing via `samtools`.
- **featureCounts (Subread) & STAR Counts**: Exon/gene-level transcriptomic quantification aggregated into a unified **Gene $\times$ Cell Expression Matrix** (`.csv`, `.tsv`, and AnnData `.h5ad`).

### 2. Downstream Single-Cell Analysis (Scanpy Ecosystem)
- **Cell & Gene Quality Control**: Automated filtering based on detected genes, total counts, and mitochondrial read percentage thresholds.
- **Normalization & Feature Selection**: Library size depth scaling, $\log(1 + x)$ variance stabilization, and Highly Variable Gene (HVG) selection.
- **Dimensionality Reduction & Graph Embedding**: PCA decomposition, neighborhood graph construction, and non-linear **UMAP** visualization.
- **Community Clustering**: **Leiden** graph-based clustering for cell state and subpopulation identification.
- **Differential Expression & Marker Discovery**: Identification of cluster-specific biomarker genes (Wilcoxon rank-sum / t-test) and export to `reports/downstream/cluster_markers.csv`.
- **Publication-Ready Visualization**: Automated generation of UMAP cluster plots, QC violin plots, and Marker DotPlots.

---

## 📁 Repository Structure

```
single_cell_pipeline/
├── pixi.toml                   # Pixi package & environment configuration
├── README.md                   # Project documentation & user guide
├── config/
│   └── pipeline_config.yaml    # Centralized pipeline & aligner/quant/downstream configuration
├── src/
│   ├── __init__.py
│   ├── utils.py                # Logging, command runner, summary tables
│   ├── sc_data_generator.py    # Synthetic Smart-seq2 multi-cell FASTQ & mini reference generator
│   ├── sc_qc.py                # FastQC and MultiQC engine
│   ├── sc_preprocess.py        # fastp paired-end trimming engine
│   ├── aligner.py              # STAR and HISAT2 splice-aware aligners + samtools
│   ├── quantifier.py           # featureCounts & STAR counts + AnnData matrix export
│   ├── downstream.py           # Scanpy downstream QC, PCA, UMAP, Leiden & Marker discovery
│   ├── pipeline.py             # Unified CLI workflow orchestration
│   └── build_notebooks.py      # Notebook generator & maintenance script
├── notebooks/
│   ├── 01_raw_qc_and_trimming.ipynb            # Raw FastQC, fastp trimming & metrics
│   ├── 02_alignment_and_quantification.ipynb   # STAR / HISAT2 mapping & featureCounts
│   ├── 03_single_cell_matrix_and_multiqc.ipynb # Count matrix QC, AnnData & MultiQC
│   └── 04_downstream_single_cell_analysis.ipynb# Scanpy QC, UMAP, Leiden & Marker discovery
├── data/
│   ├── raw/                    # Raw paired-end FASTQ inputs per cell (e.g. cell_01_R1/R2)
│   ├── clean/                  # Trimmed clean FASTQ outputs
│   ├── aligned/                # Coordinate-sorted and indexed BAM files (.bam, .bai)
│   ├── counts/                 # Gene-by-cell matrix (CSV, TSV, AnnData .h5ad)
│   └── reference/              # Reference genome FASTA, GTF, and aligner indices
└── reports/
    ├── qc_raw/                 # Raw FastQC reports
    ├── qc_clean/               # Clean FastQC reports
    ├── fastp/                  # fastp trimming HTML & JSON reports
    ├── multiqc/                # MultiQC aggregated interactive dashboard
    └── downstream/             # UMAP plots, QC violins, Marker DotPlots & CSV tables
```

---

## ⚡ Quick Start

### 1. Initialize the Environment
```bash
pixi install
```

### 2. Generate Sample Smart-seq2 Data & Mini Reference
```bash
pixi run generate-sample-data
```

### 3. Run the Full End-to-End Pipeline (Upstream + Downstream)
```bash
pixi run full-pipeline
```

---

## 🛠️ CLI Subcommands & Pixi Tasks

| Pixi Task | Description |
| :--- | :--- |
| `pixi run generate-sample-data` | Generates 8 synthetic Smart-seq2 cell FASTQ pairs + mini genome/GTF |
| `pixi run qc-raw` | Runs FastQC on raw FASTQ files |
| `pixi run trim` | Performs fastp paired-end trimming (Nextera adapter, poly-A/G) |
| `pixi run qc-clean` | Runs FastQC on clean FASTQ files |
| `pixi run align-star` | Aligns clean reads using STAR with BAM sorting & indexing |
| `pixi run align-hisat2` | Aligns clean reads using HISAT2 with BAM sorting & indexing |
| `pixi run quant-featurecounts` | Quantifies genes using Subread featureCounts |
| `pixi run quant-star` | Quantifies genes using STAR ReadsPerGene counts |
| `pixi run multiqc-report` | Aggregates all reports into an interactive MultiQC dashboard |
| `pixi run downstream` | Executes Scanpy downstream QC, PCA, UMAP, Leiden & Marker discovery |
| `pixi run full-pipeline` | Executes the complete workflow from raw FASTQs to downstream results |
| `pixi run jupyter` | Launches Jupyter Lab (`http://localhost:8888`) |
| `pixi run test-notebooks` | Programmatically executes and validates all 4 notebooks |

---

## 📓 Interactive Jupyter Notebooks Guide

Launch Jupyter Lab with:
```bash
pixi run jupyter
```

1. **`notebooks/01_raw_qc_and_trimming.ipynb`**:
   - Inspect raw paired-end FASTQ structure.
   - Run FastQC and examine per-base quality profiles.
   - Run `fastp` adapter and poly-A/G tail trimming.
   - Before vs. After read retention and Q30 score comparisons.

2. **`notebooks/02_alignment_and_quantification.ipynb`**:
   - Reference indexing and splice-aware mapping with STAR and HISAT2.
   - SAM/BAM sorting and indexing.
   - Gene quantification via `featureCounts` and STAR `GeneCounts`.
   - Unique mapping rate benchmarking.

3. **`notebooks/03_single_cell_matrix_and_multiqc.ipynb`**:
   - Load generated `gene_cell_count_matrix.h5ad` into **AnnData**.
   - Calculate single-cell QC metrics (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`).
   - Visualize cell distribution metrics.
   - Compile and view the **MultiQC** interactive report.

4. **`notebooks/04_downstream_single_cell_analysis.ipynb`**:
   - Cell & gene QC filtering and mitochondrial thresholding.
   - Depth normalization and Highly Variable Gene (HVG) selection.
   - Principal Component Analysis (PCA) & variance explained.
   - Non-linear UMAP 2D projection and Leiden cluster detection.
   - Biomarker discovery (Differential Expression) and DotPlot visualization.

---

## 🧪 Automated Testing

To run the complete verification suite including pipeline execution and automated notebook execution:
```bash
pixi run full-pipeline
pixi run test-notebooks
```
