# 🧬 Production Smart-seq2 Single-Cell Upstream Pipeline

A high-performance, reproducible bioinformatics pipeline for processing **Smart-seq2 and plate-based Single-Cell RNA-seq data** from real-world paired-end FASTQ reads to alignments, gene quantification count matrices, and interactive **MultiQC** quality dashboards.

> **Note**: This repository focuses strictly on **production upstream processing** (Raw FASTQ $\rightarrow$ Count Matrix & MultiQC). All downstream tasks (Scanpy, UMAP, clustering) are excluded.

---

## 🏛️ Dual-Track Architecture

This project is divided into two completely independent workflows:

```
single_cell_pipeline/
├── PLAN.md                     # Master AI roadmap, context strategy & milestones
├── pixi.toml                   # Unified environment (BioConda / Conda-Forge)
├── config.env.example          # Production configuration template
├── data/
│   ├── raw/                    # User-supplied raw paired-end FASTQ reads
│   ├── reference/              # Reference genome FASTA, GTF, and STAR index
│   ├── clean/                  # Trimmed FASTQ reads
│   ├── aligned/                # Coordinate-sorted BAMs & indices (.bam, .bai)
│   └── counts/                 # Gene-by-cell count matrix (.tsv)
│
├── scripts/                    # [Part 2.1] Pure Bash Modular Tool Suite
│   ├── env.sh                  # Common environment variables & hardware detection
│   ├── setup_reference.sh      # Automated reference download & STAR indexing
│   ├── qc_data_00.sh           # Stage 00: FastQC raw reads
│   ├── fastp_01.sh             # Stage 01: fastp adapter & quality trimming
│   ├── align_02.sh             # Stage 02: STAR alignment & samtools indexing
│   ├── quant_03.sh             # Stage 03: featureCounts & AWK matrix assembly
│   ├── multiqc_04.sh           # Stage 04: MultiQC report aggregation
│   └── run_all.sh              # Master runner with CLI flags & validation
│
└── nextflow/                   # [Part 2.2] Production Nextflow DSL2 Pipeline
    ├── main.nf                 # Top-level workflow orchestration
    ├── nextflow.config         # Profiles (local, docker, singularity, slurm), BioContainers
    ├── samplesheet.example.csv # Template samplesheet CSV
    └── modules/                # Self-contained DSL2 tool modules
        ├── fastqc.nf           # FastQC process
        ├── fastp.nf            # fastp trimming process
        ├── star.nf             # STAR alignment & BAM indexing process
        ├── featurecounts.nf    # featureCounts quantification process
        └── multiqc.nf          # MultiQC aggregation process
```

---

## 🚀 Getting Started

### 1. Environment Setup via [Pixi](https://pixi.sh/)
```bash
pixi install
```

### 2. Set Up Reference Genome (Human or Mouse)
Use the automated reference downloader and indexer:
```bash
# Download GENCODE Human (GRCh38) and build STAR index:
bash scripts/setup_reference.sh --species human --threads 16

# Or download GENCODE Mouse (GRCm39):
bash scripts/setup_reference.sh --species mouse --threads 16

# Or index custom local genome files:
bash scripts/setup_reference.sh --fasta /path/to/genome.fa --gtf /path/to/genes.gtf --threads 16
```

---

## 💻 Track 1: Pure Bash Modular Scripts (Part 2.1)

Designed for **maximum parameter transparency**, tool inspection, and manual execution. **Zero Python dependencies are used**.

### Run with CLI arguments:
```bash
bash scripts/run_all.sh \
    --indir /path/to/fastqs \
    --star-index data/reference/star_index \
    --gtf data/reference/genes.gtf \
    --threads 16
```

### Or run individual stages:
```bash
bash scripts/qc_data_00.sh data/raw reports/qc_raw 16
bash scripts/fastp_01.sh data/raw data/clean reports/fastp 16
bash scripts/align_02.sh data/clean data/aligned data/reference/star_index 16
bash scripts/quant_03.sh data/aligned data/counts data/reference/genes.gtf 16
bash scripts/multiqc_04.sh reports/multiqc
```

---

## 🌊 Track 2: Production Nextflow DSL2 Pipeline (Part 2.2)

Designed for **enterprise deployment**, massive parallelism, cluster execution (Slurm), and cloud portability via BioContainers. Does **not** depend on Track 1 scripts.

### Run via Samplesheet:
```bash
# Create samplesheet.csv:
# sample,fastq_1,fastq_2
# cell_01,/data/cell_01_R1.fastq.gz,/data/cell_01_R2.fastq.gz

nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --star_index data/reference/star_index \
    --gtf data/reference/genes.gtf \
    --outdir results \
    -profile docker
```

### Run on HPC with Singularity:
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --star_index data/reference/star_index \
    --gtf data/reference/genes.gtf \
    --outdir results \
    -profile slurm,singularity
```

---

## 📋 Roadmap & Milestone Plan

For detailed context management rules, tool arguments cheatsheets, and AI guidelines, see **[`PLAN.md`](PLAN.md)**.
