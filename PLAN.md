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

---

## 📖 Tool Argument Anatomy & Cheatsheet

### 1. `fastp` (Trimming)
```bash
fastp \
  --in1 sample_R1_001.fastq.gz --in2 sample_R2_001.fastq.gz \
  --out1 clean_R1.fastq.gz --out2 clean_R2.fastq.gz \
  --detect_adapter_for_pe \                      # Auto-detect paired adapters
  --trim_poly_g --poly_g_min_len 10 \            # Remove NovaSeq poly-G artifacts
  --trim_poly_x --poly_x_min_len 10 \            # Remove poly-A/T/C homopolymers
  --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
  --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
  -q 20 -u 30 \                                  # Min Q20, max 30% unqualified bases
  -l 35 \                                        # Production length cutoff (>= 35 bp)
  --thread 8 --json report.json --html report.html
```

### 2. `STAR` (Mammalian Alignment)
```bash
STAR \
  --runThreadN 16 \
  --genomeDir /refs/GRCh38_star_index \
  --readFilesIn clean_R1.fastq.gz clean_R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM NM MD jM jI XS \
  --outFilterType BySJout \
  --limitBAMsortRAM 31000000000 \                # Prevent OOM during sorting
  --quantMode GeneCounts \
  --outFileNamePrefix data/aligned/sample_
```

### 3. `subread featureCounts` (Quantification)
```bash
featureCounts \
  -T 16 \
  -p --countReadPairs \                          # Paired-end fragment counting
  -t exon -g gene_id \                           # Gene-level summarization
  -s 0 \                                         # 0 = unstranded (Smart-seq2)
  -Q 10 \                                        # Min MAPQ 10
  -a /refs/gencode.v44.gtf \
  -o data/counts/featurecounts_raw.txt \
  data/aligned/*.sortedByCoord.out.bam
```

---

## 🤝 Handoff Protocol for AI Sessions

When concluding an AI turn or handing off to another agent, always provide:
1. **Target Milestone**: Which milestone from Section 6 was worked on.
2. **Commands Executed**: Exact bash/nextflow commands tested.
3. **Artifacts Produced**: Output files created.
4. **Next Immediate Steps**: Concrete single task for the next session.
