# 🧬 Single-Cell RAW to Clean FASTQ Preprocessing & QC Suite

A production-grade, reproducible pipeline and interactive Jupyter analysis environment for processing **Single-Cell RNA-seq FASTQ** data (e.g. 10x Genomics Chromium 3' v2/v3, Drop-seq, Split-seq) from raw sequencing reads to quality-controlled, adapter-trimmed, error-corrected clean FASTQ files.

Managed via **[Pixi](https://pixi.sh/)** with native `bioconda` and `conda-forge` environments.

---

## 🚀 Key Features

- **⚡ Complete Environment Reproducibility**: Zero conda-activation headaches; managed entirely through `pixi.toml` with fast C++/Rust resolution.
- **🔬 Comprehensive Single-Cell QC**:
  - Per-cycle Phred quality profiles for Cell Barcode, UMI, and cDNA.
  - Cell Barcode rank distribution and **Knee Plot** cell vs. ambient droplet inflection estimation.
  - Barcode Whitelist match rates (Exact match vs. 1-Hamming distance recoverable vs. invalid).
  - FastQC and automated **MultiQC** interactive HTML reporting.
- **✂️ Robust Preprocessing & Trimming**:
  - Automated adapter / Template Switch Oligo (TSO) removal.
  - Poly-G (NextSeq/NovaSeq two-color artifact) and Poly-A tail trimming via `fastp`.
  - Sliding-window low-quality 3' tail trimming.
- **🎯 1-Hamming Distance Barcode Error Correction**:
  - Recovers sequencing errors in cell barcodes against official 10x whitelists.
  - Discards uncorrectable chimeric and low-quality barcodes.
- **📓 3 Interactive Jupyter Notebooks**:
  - Step-by-step visual exploration, parameter tuning, and publication-ready plots.
- **📦 Downstream Compatibility**:
  - Produces clean paired FASTQs (`*_val_R1.fastq.gz`, `*_val_R2.fastq.gz`) and header-extracted FASTQs (`@READ_CB_UMI`) compatible with **STARsolo**, **Kallisto/Bustools**, **CellRanger**, and **Salmon/Alevin**.

---

## 📁 Repository Structure

```
single_cell_pipeline/
├── pixi.toml                   # Pixi package & environment configuration
├── README.md                   # Project documentation & user guide
├── config/
│   └── pipeline_config.yaml    # Centralized pipeline & chemistry configuration
├── src/
│   ├── __init__.py
│   ├── utils.py                # Logging, command runner, summary tables
│   ├── sc_data_generator.py    # Synthetic 10x FASTQ & whitelist generator
│   ├── sc_qc.py                # Single-cell QC & Knee plot analysis engine
│   ├── sc_preprocess.py        # fastp trimming & 1-Hamming barcode error corrector
│   ├── pipeline.py             # Unified CLI workflow orchestration
│   └── build_notebooks.py      # Notebook generator & maintenance script
├── notebooks/
│   ├── 01_raw_data_qc.ipynb                # Raw QC, per-cycle Phred, Knee plot
│   ├── 02_preprocessing_pipeline.ipynb     # Interactive trimming & error-correction
│   └── 03_post_qc_and_benchmarking.ipynb   # Before vs After benchmarking & MultiQC
├── data/
│   ├── raw/                    # Raw paired FASTQ inputs
│   ├── clean/                  # Preprocessed clean FASTQ outputs
│   └── whitelist/              # 10x barcode whitelists (e.g. 737K-august-2016.txt)
└── reports/
    ├── qc_raw/                 # Raw FastQC & JSON metrics
    ├── qc_clean/               # Clean FastQC & JSON metrics
    ├── fastp/                  # fastp trimming HTML & JSON reports
    └── multiqc/                # MultiQC aggregated interactive dashboard
```

---

## ⚡ Quick Start

### 1. Install Pixi (if not already installed)
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Initialize the Environment
```bash
pixi install
```

### 3. Generate Sample 10x scRNA-seq Data (or use your own)
```bash
pixi run generate-sample-data
```
*Generates realistic 10x Chromium 3' v3 paired FASTQs with cell barcodes, UMIs, adapter contamination, low-quality ends, poly-A/poly-G tails, and an authentic 10x whitelist.*

### 4. Run the Full End-to-End Pipeline
```bash
pixi run full-pipeline
```

---

## 🛠️ CLI Subcommands & Pixi Tasks

| Pixi Task | Description |
| :--- | :--- |
| `pixi run generate-sample-data` | Simulates realistic 10x single-cell FASTQs |
| `pixi run qc-raw` | Computes raw FastQC & single-cell QC metrics |
| `pixi run preprocess` | Performs fastp trimming + whitelist error correction |
| `pixi run qc-clean` | Computes clean FastQC & single-cell QC metrics |
| `pixi run multiqc-report` | Aggregates all reports into an interactive MultiQC dashboard |
| `pixi run full-pipeline` | Runs the complete workflow from start to finish |
| `pixi run jupyter` | Launches Jupyter Lab (`http://localhost:8888`) |
| `pixi run test-notebooks` | Programmatically executes and validates all 3 notebooks |

You can also run individual steps directly via the CLI:
```bash
# Run QC on raw FASTQs
pixi run python -m src.pipeline qc --stage raw --config config/pipeline_config.yaml

# Run Preprocessing
pixi run python -m src.pipeline preprocess --config config/pipeline_config.yaml

# Run QC on clean FASTQs
pixi run python -m src.pipeline qc --stage clean --config config/pipeline_config.yaml

# Build MultiQC report
pixi run python -m src.pipeline multiqc --config config/pipeline_config.yaml
```

---

## 📓 Interactive Jupyter Notebooks Guide

Launch Jupyter with:
```bash
pixi run jupyter
```

1. **`notebooks/01_raw_data_qc.ipynb`**:
   - Inspect raw FASTQ structure (R1 28bp CB+UMI, R2 91bp cDNA).
   - Per-cycle Phred scores across barcode cycles.
   - Barcode frequency distribution and interactive **Knee Plot** with automatic cell inflection cutoff.
   - Whitelist match breakdown (Exact vs 1-bp mismatch vs invalid).

2. **`notebooks/02_preprocessing_pipeline.ipynb`**:
   - Interactive `fastp` adapter, poly-G, and poly-A trimming.
   - 1-Hamming distance error correction of mutated barcodes.
   - Yield breakdown and clean FASTQ generation.

3. **`notebooks/03_post_qc_and_benchmarking.ipynb`**:
   - Side-by-side Before vs. After comparison (Q30 gain, read loss reasons).
   - Knee Plot before vs. after error correction.
   - Embedded interactive **MultiQC** dashboard.
   - Downstream alignment commands for STARsolo, Kallisto, and CellRanger.

---

## ⚙️ Custom Configuration (`config/pipeline_config.yaml`)

Edit `config/pipeline_config.yaml` to customize chemistry, trimming thresholds, and file paths:

```yaml
# Library Chemistry
chemistry:
  name: "10x_v3"
  r1_structure:
    cell_barcode_len: 16
    umi_len: 12
    total_len: 28

# Preprocessing & Trimming Parameters
preprocessing:
  threads: 4
  fastp:
    qualified_quality_phred: 20
    unqualified_percent_limit: 30
    min_length: 25
    trim_poly_g: true
    trim_poly_x: true
    cut_front: true
    cut_tail: true
    adapter_sequence_r2: "AAGCAGTGGTATCAACGCAGAGTACATGGG" # 10x TSO
  
  barcode_filtering:
    max_hamming_distance: 1
    filter_unknown_barcodes: true
```

---

## 🎯 Downstream Alignment Examples

Once cleaned, you can feed the output files directly into standard alignment tools:

### STARsolo
```bash
STAR --genomeDir /path/to/star_index \
     --readFilesIn data/clean/scRNA_sample_01_val_R2.fastq.gz data/clean/scRNA_sample_01_val_R1.fastq.gz \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist data/whitelist/737K-august-2016.txt \
     --soloCBlen 16 --soloUMIlen 12 \
     --soloStrand Forward \
     --outSAMtype BAM SortedByCoordinate
```

### Kallisto / Bustools
```bash
kallisto bus -i transcripts.idx -o bus_output/ -x 10xv3 -t 4 \
  data/clean/scRNA_sample_01_val_R1.fastq.gz data/clean/scRNA_sample_01_val_R2.fastq.gz
bustools sort -t 4 -o bus_output/sorted.bus bus_output/output.bus
bustools count -o bus_output/count -g transcripts_to_genes.txt -e bus_output/matrix.ec -t bus_output/transcripts.txt --genecounts bus_output/sorted.bus
```

---

## 🧪 Automated Testing

To run the complete verification suite including pipeline execution and automated notebook execution:
```bash
pixi run full-pipeline
pixi run test-notebooks
```
