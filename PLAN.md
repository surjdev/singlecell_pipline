# 🧬 Production Smart-seq2 Single-Cell Upstream Pipeline: Roadmap & Milestone Plan (`PLAN.md`)

> **Project Goal**: Production-grade upstream processing of Smart-seq2 and plate-based single-cell RNA-seq data from real-world paired-end FASTQ reads to gene-by-cell count matrices and MultiQC quality dashboards. All downstream tasks (Scanpy, clustering, UMAP, marker discovery) are strictly excluded from this repository.

---

## 📑 Table of Contents
1. [AI Context Management Protocol](#-ai-context-management-protocol-critical-for-future-ai-turns)
2. [Dual-Track Architecture Overview](#-dual-track-architecture-overview)
3. [Production Reference Preparation](#-production-reference-preparation)
4. [Part 2.1: Pure Bash Modular Tool Suite (Production Grade)](#-part-21-pure-bash-modular-tool-suite)
5. [Part 2.2: Production Nextflow DSL2 Pipeline (Enterprise Ready)](#-part-22-production-nextflow-dsl2-pipeline)
6. [Milestones & Progress Checklist](#-milestones--progress-checklist)
7. [Tool Argument Anatomy & Cheatsheet](#-tool-argument-anatomy--cheatsheet)
8. [Handoff Protocol for AI Sessions](#-handoff-protocol-for-ai-sessions)

---

## 🧠 AI Context Management Protocol (Critical for Future AI Turns)

To prevent context exhaustion and maintain high coding accuracy across iterative AI sessions, all subsequent AI agents **MUST** adhere to the following rules:

### 1. Token Economy & Scoping Rules
- **One Component per Turn**: Focus exclusively on one module, script, or configuration at a time. Never load the entire codebase into context simultaneously.
- **Do NOT Dump Large Output Files**:
  - Never dump full FASTQ, SAM, BAM, GTF, or raw count matrix files into chat context.
  - Inspect files using bounded commands:
    ```bash
    head -n 10 data/counts/gene_cell_count_matrix.tsv
    samtools view -h data/aligned/sample_01_Aligned.sortedByCoord.out.bam | head -n 20
    samtools flagstat data/aligned/sample_01_Aligned.sortedByCoord.out.bam
    ```
- **Check Status Before Acting**: Run `git status` and inspect `PLAN.md` milestone status before initiating changes.

### 2. Strict Input/Output Contracts Between Modules

| Stage | Input Contract | Output Contract | Log / Report Contract |
|---|---|---|---|
| **00. Raw QC** | `data/raw/*_R{1,2}*.fastq.gz` | None (QC only) | `reports/qc_raw/*.{html,zip}` |
| **01. Trimming** | `data/raw/*_R{1,2}*.fastq.gz` | `data/clean/*_R{1,2}.clean.fastq.gz` | `reports/fastp/*.{json,html}` |
| **02. Alignment** | `data/clean/*_R{1,2}.clean.fastq.gz` + `STAR_INDEX` | `data/aligned/*_Aligned.sortedByCoord.out.bam` + `.bai` | `data/aligned/*_Log.final.out` |
| **03. Quant** | `data/aligned/*_Aligned.sortedByCoord.out.bam` + `GTF_FILE` | `data/counts/gene_cell_count_matrix.tsv` | `data/counts/featurecounts_raw.txt.summary` |
| **04. MultiQC** | All report directories above | None | `reports/multiqc/single_cell_multiqc_report.html` |

---

## 🏗️ Dual-Track Architecture Overview

```
single_cell_pipeline/
├── PLAN.md                     # Master AI roadmap & context strategy (this file)
├── pixi.toml                   # Unified environment (BioConda/Conda-forge)
├── config.env.example          # Production configuration template
├── data/
│   ├── raw/                    # User-supplied raw paired-end FASTQ reads
│   ├── reference/              # Reference genome FASTA, GTF, and STAR index
│   ├── clean/                  # Trimmed FASTQ reads
│   ├── aligned/                # Coordinate-sorted BAMs & indices
│   └── counts/                 # Gene-by-cell count matrix (TSV)
│
├── scripts/                    # [Part 2.1] Pure Bash Modular Tool Suite
│   ├── env.sh                  # Common environment variables & hardware detection
│   ├── setup_reference.sh      # Automated reference download & STAR indexing
│   ├── qc_data_00.sh           # Stage 00: FastQC on raw reads
│   ├── fastp_01.sh             # Stage 01: fastp adapter & quality trimming
│   ├── align_02.sh             # Stage 02: STAR alignment & samtools indexing
│   ├── quant_03.sh             # Stage 03: featureCounts & AWK matrix formatting
│   ├── multiqc_04.sh           # Stage 04: MultiQC log aggregation
│   └── run_all.sh              # Master runner with CLI flags & validation
│
└── nextflow/                   # [Part 2.2] Production Nextflow DSL2 Pipeline
    ├── main.nf                 # Top-level workflow orchestration (samplesheet / glob)
    ├── nextflow.config         # Configuration, profiles, resource labels & BioContainers
    ├── samplesheet.example.csv # Template samplesheet CSV
    └── modules/                # Self-contained DSL2 tool modules
        ├── fastqc.nf           # FastQC process
        ├── fastp.nf            # fastp trimming process
        ├── star.nf             # STAR alignment process
        ├── featurecounts.nf    # featureCounts quantification process
        └── multiqc.nf          # MultiQC aggregation process
```

> [!IMPORTANT]
> **Track Independence**:
> - **Part 2.1 (Bash)**: Provides maximum transparency and direct parameter inspection with pure Bash without Python.
> - **Part 2.2 (Nextflow)**: Provides enterprise-level scalability, cloud/cluster portability via BioContainers, and automatic failure recovery. Does **not** call Part 2.1 scripts.

---

## 🧬 Production Reference Preparation

To prepare a production reference genome and STAR index (e.g. human GRCh38 or mouse GRCm39):

```bash
# Automated download (GENCODE) and STAR indexing for Human:
bash scripts/setup_reference.sh --species human --threads 16

# Automated download (GENCODE) and STAR indexing for Mouse:
bash scripts/setup_reference.sh --species mouse --threads 16

# Or build index from existing local FASTA and GTF:
bash scripts/setup_reference.sh \
    --fasta /path/to/GRCh38.genome.fa \
    --gtf /path/to/gencode.v44.gtf \
    --threads 16
```

---

## 🔧 Part 2.1: Pure Bash Modular Tool Suite

Each script in `scripts/` is completely standalone, accepts environment overrides, uses `set -euo pipefail` for strict error handling, and formats outputs using native UNIX utilities (`awk`, `cut`, `sed`).

### Running Real Data:
```bash
# 1. (Optional) Customize configuration:
cp config.env.example config.env
nano config.env

# 2. Run the end-to-end production bash pipeline:
bash scripts/run_all.sh \
    --indir /path/to/real_fastqs \
    --star-index /path/to/star_index \
    --gtf /path/to/genes.gtf \
    --threads 16

# Or execute individual stages:
bash scripts/qc_data_00.sh /path/to/real_fastqs reports/qc_raw 16
bash scripts/fastp_01.sh /path/to/real_fastqs data/clean reports/fastp 16
bash scripts/align_02.sh data/clean data/aligned /path/to/star_index 16
bash scripts/quant_03.sh data/aligned data/counts /path/to/genes.gtf 16
bash scripts/multiqc_04.sh reports/multiqc
```

---

## 🚀 Part 2.2: Production Nextflow DSL2 Pipeline

Built with standard Nextflow DSL2 syntax, channel operators, dynamic resource scaling, and BioContainer support.

### Running Real Data with Samplesheet:
```bash
# Prepare samplesheet.csv:
# sample,fastq_1,fastq_2
# cell_01,/data/cell_01_R1.fastq.gz,/data/cell_01_R2.fastq.gz

# Execute via Docker:
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --star_index /path/to/star_index \
    --gtf /path/to/genes.gtf \
    --outdir results_production \
    -profile docker

# Execute on HPC / Slurm cluster with Singularity:
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --star_index /path/to/star_index \
    --gtf /path/to/genes.gtf \
    --outdir results_production \
    -profile slurm,singularity
```

---

## 📊 Milestones & Progress Checklist

### Milestone 0: Clean Up & Test Elimination
- [x] Delete all legacy downstream Python code and notebooks.
- [x] Delete all toy/synthetic test data, mock references, and toy caches from workspace.
- [x] Create clean `.gitkeep` directory skeleton (`data/raw`, `data/reference`, etc.).

### Milestone 1: Production Bash Tool Suite (Part 2.1)
- [x] Configurable `scripts/env.sh` with hardware auto-detection and RAM bounds.
- [x] `scripts/setup_reference.sh`: Automated GENCODE human/mouse download and STAR indexer.
- [x] `scripts/qc_data_00.sh`: Batch FastQC with xargs support.
- [x] `scripts/fastp_01.sh`: Multi-pattern Illumina pair matching, poly-G/X, and configurable read length.
- [x] `scripts/align_02.sh`: Production STAR alignment with memory safety limits (`--limitBAMsortRAM`).
- [x] `scripts/quant_03.sh`: Subread featureCounts with configurable strandedness and pure AWK matrix formatter.
- [x] `scripts/multiqc_04.sh`: Multi-source log aggregation.
- [x] `scripts/run_all.sh`: Production driver with CLI argument parsing.

### Milestone 2: Production Nextflow DSL2 Pipeline (Part 2.2)
- [x] `nextflow/nextflow.config`: Production resource management, BioContainers, and profiles.
- [x] Samplesheet CSV parsing support (`sample,fastq_1,fastq_2`).
- [x] `nextflow/main.nf`: Robust reference validation and DSL2 dataflow.
- [x] BioContainers integrated for FastQC, fastp, STAR, Subread, MultiQC.

### Milestone 3: Downstream Analysis Workbench (`workbench/` Sub-Project)
- [x] **Decoupled Architecture**: Completely self-contained sub-project with its own `pyproject.toml`, `pixi.toml`, `.gitignore`, and `README.md` (ready for independent repo extraction).
- [x] **Strict Classical Focus (No DL)**: Deep learning modules explicitly rolled back per user command ("Rollback ไม่เอา DL แล้ว"). Zero PyTorch / scvi-tools overhead.
- [x] **`sc_workbench` Library**:
  - `io.py`: Ingestion of upstream `data/counts/gene_cell_count_matrix.tsv`, 10x MTX, H5AD, and Seurat CSV export.
  - `qc.py`: Total counts, detected genes, mito/ribo %, Scrublet doublet detection, filtering.
  - `preprocess.py`: Library size normalization, log1p, HVG selection, scaling.
  - `reduction.py`: PCA, Harmony multi-batch integration, k-NN graph, UMAP, t-SNE.
  - `clustering.py`: Leiden (igraph) & Louvain community detection, sub-clustering.
  - `markers.py`: Wilcoxon / Welch differential expression, tidy DataFrame extraction.
  - `annotation.py`: Gene signature scoring and cell type annotation.
  - `pathway.py`: GSEAPY Enrichr over-representation analysis & GSEA.
  - `trajectory.py`: PAGA connectivity graph and Diffusion Pseudotime (DPT).
  - `plotting.py`: Publication palettes, QC violins, UMAP, DotPlot, Volcano plots.
  - `workbench.py`: High-level `SingleCellWorkbench` fluent chaining pipeline.
- [x] **Jupyter Notebook Suite (`workbench/notebooks/`)**:
  - `01_ingest_and_qc.ipynb`: Upstream TSV ingest, QC metrics, Scrublet doublet filtering.
  - `02_clustering_and_umap.ipynb`: Normalization, HVG, PCA, Harmony, kNN, UMAP, Leiden.
  - `03_marker_genes_and_annotation.ipynb`: Differential expression, DotPlots, cell annotation.
  - `04_pathway_and_trajectory.ipynb`: Enrichr pathways, PAGA lineage, DPT pseudotime.
- [x] **Automated Testing Suite (`workbench/tests/`)**:
  - `test_io.py`, `test_qc.py`, `test_pipeline.py` (100% passing in Pixi).

---

## 📖 Tool Argument Anatomy, Error Troubleshooting & Output Interpretation

> An exhaustive operational manual for every tool is available in **[`docs/TOOLS_RUNBOOK.md`](docs/TOOLS_RUNBOOK.md)**. Below is a structured summary covering all 6 upstream tools:

### 1. FastQC (Quality Assessment)
- **Run Command**:
  ```bash
  fastqc --outdir reports/qc_raw --threads 8 --noextract --quiet data/raw/*.fastq.gz
  ```
- **Flags Breakdown**:
  - `--outdir / -o`: Output directory for HTML and ZIP results.
  - `--threads / -t`: Number of parallel files to analyze simultaneously (~250MB RAM/thread).
  - `--noextract`: Prevents uncompressing ZIP file, saving disk storage.
  - `--quiet / -q`: Suppresses verbose stdout messages during pipeline runs.
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `java.lang.OutOfMemoryError: Java heap space`: Exceeded default 512MB heap. Fix: `export _JAVA_OPTIONS="-Xmx2048m"`.
  - `gzip: unexpected end of file`: Corrupted/truncated FASTQ file. Fix: Verify with `gzip -t <file>` and re-download.
  - `Too many open files`: Ulimit exceeded. Fix: Batch with `xargs -n 20`.
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `<sample>_fastqc.html`: Interactive visual report.
  - `<sample>_fastqc.zip`: Raw metric data table (`fastqc_data.txt`).
  - *Healthy standard*: Phred quality score > 28, GC content unimodal bell curve, adapter content drops to 0% after fastp.

---

### 2. fastp (Trimming & Quality Filtering)
- **Run Command**:
  ```bash
  fastp \
    --in1 data/raw/sample_R1.fastq.gz --in2 data/raw/sample_R2.fastq.gz \
    --out1 data/clean/sample_R1.clean.fastq.gz --out2 data/clean/sample_R2.clean.fastq.gz \
    --detect_adapter_for_pe \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x --poly_x_min_len 10 \
    --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
    --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
    -q 20 -u 30 -l 35 --thread 8 \
    --json reports/fastp/sample_fastp.json --html reports/fastp/sample_fastp.html
  ```
- **Flags Breakdown**:
  - `--detect_adapter_for_pe`: Automatically detects paired adapter sequences via read overlap.
  - `--adapter_sequence / --adapter_sequence_r2`: Custom sequence override (e.g. Nextera `CTGTCTCTTATACACATCT`).
  - `--trim_poly_g`: Trims NovaSeq/NextSeq 2-color dark cycle poly-G artifacts.
  - `--trim_poly_x`: Trims cDNA poly-A/T/C homopolymers.
  - `--cut_front / --cut_tail`: 5' and 3' sliding window trimming based on mean Phred quality.
  - `-q 20 -u 30`: Base quality threshold (Q20) and max allowed unqualified base percentage (30%).
  - `-l 35`: Minimum read length cutoff (filters out short ambiguous reads).
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `R1 and R2 have different number of reads`: Pairing mismatch or truncated mate. Fix: check line count symmetry with `zcat <file> | wc -l`.
  - `> 90% reads filtered out`: `-l` length requirement set higher than actual sequencer cycle length. Fix: check raw length in FastQC and lower `-l`.
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `sample_R{1,2}.clean.fastq.gz`: Cleaned paired-end reads.
  - `sample_fastp.json` & `.html`: Trimming statistics and before/after quality curves.
  - *Healthy standard*: > 85-95% reads passed filter, post-trim Q30 > 90-95%.

---

### 3. STAR (Splice-Aware Alignment)
- **Run Command**:
  ```bash
  STAR \
    --runThreadN 16 \
    --genomeDir data/reference/star_index \
    --readFilesIn data/clean/sample_R1.clean.fastq.gz data/clean/sample_R2.clean.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM NM MD jM jI XS \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --limitBAMsortRAM 31000000000 \
    --quantMode GeneCounts \
    --outFileNamePrefix data/aligned/sample_
  ```
- **Flags Breakdown**:
  - `--genomeDir`: Directory containing pre-built STAR suffix array index.
  - `--readFilesCommand zcat`: Decompresses gzipped input FASTQ on the fly.
  - `--outSAMtype BAM SortedByCoordinate`: Directly produces coordinate-sorted BAM (no extra samtools sort step).
  - `--outSAMunmapped Within`: Retains unmapped reads inside BAM file to preserve total library read counts.
  - `--limitBAMsortRAM 31000000000`: Buffer size for coordinate sorting (~31GB). Prevents OOM crashes on mammalian genomes.
  - `--outFilterMultimapNmax 20`: Maximum alignment loci allowed before discarding multimappers.
  - `--quantMode GeneCounts`: STAR internal read summarization per gene.
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `BAMsort: limitBAMsortRAM reached / Killed / Exit 137`: Out of memory during BAM sorting. Fix: increase `--limitBAMsortRAM` or decrease `--runThreadN`.
  - `FATAL INPUT ERROR: could not open genome file ... Genome`: Corrupted or missing STAR index. Fix: build with `scripts/setup_reference.sh`.
  - `Uniquely mapped reads % < 40%`: Wrong species reference or severe contamination. Fix: verify sample organism.
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `<sample>_Aligned.sortedByCoord.out.bam`: Coordinate-sorted BAM file.
  - `<sample>_Log.final.out`: Alignment report summary.
  - `<sample>_ReadsPerGene.out.tab`: Raw gene counts table.
  - *Healthy standard*: Uniquely mapped reads > 70-85%, mismatch rate < 0.5-1%.

---

### 4. samtools (BAM Indexing & Inspection)
- **Run Command**:
  ```bash
  samtools index -@ 8 data/aligned/sample_Aligned.sortedByCoord.out.bam
  samtools flagstat data/aligned/sample_Aligned.sortedByCoord.out.bam
  ```
- **Flags Breakdown**:
  - `index -@ <threads>`: Parallel index builder producing `.bai` index.
  - `flagstat`: Alignment summary statistics.
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `[E::hts_idx_push] NO_COOR reads not in a single block at the end`: Attempting to index an unsorted or queryname-sorted BAM. Fix: sort by coordinate first.
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `<sample>_Aligned.sortedByCoord.out.bam.bai`: Coordinate index enabling fast random querying in featureCounts and IGV.

---

### 5. Subread featureCounts (Gene Quantification)
- **Run Command**:
  ```bash
  featureCounts \
    -T 16 -p --countReadPairs \
    -t exon -g gene_id \
    -s 0 -Q 10 \
    -a data/reference/genes.gtf \
    -o data/counts/featurecounts_raw.txt \
    data/aligned/*.sortedByCoord.out.bam
  ```
- **Flags Breakdown**:
  - `-p --countReadPairs`: Paired-end mode; counts fragments (pairs) rather than individual reads twice.
  - `-t exon`: Feature type to quantify in GTF column 3.
  - `-g gene_id`: Meta-feature grouping attribute in GTF column 9.
  - `-s 0`: Strand specificity (0 = unstranded Smart-seq2, 1 = stranded, 2 = reverse stranded).
  - `-Q 10`: Minimum MAPQ mapping quality threshold.
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `Assigned: 0 (0.0%)`: Chromosome naming mismatch between BAM (`chr1, chr2`) and GTF (`1, 2`) or inverted strand (`-s 1` instead of `-s 0`). Fix: ensure consistent chromosome prefixes.
  - `Failed to open annotation file`: GTF is gzipped (`.gtf.gz`). Fix: decompress GTF (`gzip -d genes.gtf.gz`).
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `featurecounts_raw.txt.summary`: Crucial breakdown of Assigned vs Unassigned reads.
  - `gene_cell_count_matrix.tsv`: Clean Gene $\times$ Cell count matrix.
  - *Healthy standard*: Assigned reads > 60-80%.

---

### 6. MultiQC (Log Aggregation & Dashboard)
- **Run Command**:
  ```bash
  multiqc --outdir reports/multiqc --filename single_cell_multiqc_report.html --force --interactive \
    reports/qc_raw reports/fastp data/aligned data/counts
  ```
- **Flags Breakdown**:
  - `--force / -f`: Overwrite existing report.
  - `--interactive`: Enables interactive web plots for zooming and filtering.
- **Common Errors & Diagnostics ("เจอ error อะไร")**:
  - `No analysis results found`: Incorrect search paths or non-standard file names. Fix: specify folders directly.
- **Expected Outputs & Interpretation ("ผลลัพธ์เป็นอย่างไร")**:
  - `single_cell_multiqc_report.html`: Comprehensive dashboard consolidating FastQC, fastp, STAR, and featureCounts metrics across all cells.

---

## 🤝 Handoff Protocol for AI Sessions

When concluding an AI turn or handing off to another agent, always provide:
1. **Target Milestone**: Which milestone from Section 6 was worked on.
2. **Commands Executed**: Exact bash/nextflow commands tested.
3. **Artifacts Produced**: Output files created.
4. **Next Immediate Steps**: Concrete single task for the next session.

