# 📖 Production Single-Cell Pipeline: Complete Tool Runbook & Troubleshooting Guide

This runbook provides an exhaustive reference for every bioinformatics tool in the upstream Smart-seq2 single-cell pipeline. For each tool, you will find:
1. **How to run it** (CLI syntax and production execution)
2. **Complete flag breakdown** (all parameters explained)
3. **Error diagnostics & troubleshooting** (common errors, causes, and concrete fixes)
4. **Expected outputs & metric interpretation** (file formats and quality benchmarks)

---

## 📑 Table of Contents
- [1. FastQC (Raw & Clean Read Quality Control)](#1-fastqc-raw--clean-read-quality-control)
- [2. fastp (Adapter, Poly-G/X, and Quality Trimming)](#2-fastp-adapter-poly-gx-and-quality-trimming)
- [3. STAR (Splice-Aware Reference Alignment)](#3-star-splice-aware-reference-alignment)
- [4. samtools (BAM Sorting & Coordinate Indexing)](#4-samtools-bam-sorting--coordinate-indexing)
- [5. Subread featureCounts (Gene Quantification)](#5-subread-featurecounts-gene-quantification)
- [6. MultiQC (Comprehensive Metric Aggregation)](#6-multiqc-comprehensive-metric-aggregation)

---

## 1. FastQC (Raw & Clean Read Quality Control)

### 📌 1.1 Purpose
FastQC provides quality control checks on raw sequence data coming from high-throughput sequencing pipelines. It evaluates per-base sequence quality, adapter contamination, GC bias, and duplication rates.

### 💻 1.2 How to Run
```bash
# Standalone execution on paired-end reads:
fastqc \
    --outdir reports/qc_raw \
    --threads 8 \
    --noextract \
    --quiet \
    data/raw/sample_01_R1.fastq.gz data/raw/sample_01_R2.fastq.gz

# Batch execution using xargs (production scale):
find data/raw -name "*.fastq.gz" | xargs -n 20 -P 1 fastqc \
    --outdir reports/qc_raw \
    --threads 8 \
    --noextract \
    --quiet
```

### 🏷️ 1.3 Flag Breakdown
| Flag | Alternative | Description | Recommended Value | Impact |
|---|---|---|---|---|
| `--outdir` | `-o` | Output directory for HTML and ZIP reports | `reports/qc_raw` | Must exist or be created beforehand |
| `--threads` | `-t` | Number of files processed simultaneously | `4 - 16` | FastQC allocates ~250MB RAM per thread |
| `--noextract` | | Do not uncompress the output ZIP archive | Enabled | Saves disk inodes and storage |
| `--quiet` | `-q` | Suppress progress messages | Enabled in pipelines | Keeps batch logs clean |
| `--format` | `-f` | Force input file format | Auto-detected | Supports `fastq`, `bam`, `sam` |
| `--contaminants` | `-c` | Custom file containing contaminant sequences | Optional | Identifies vector/primer sequences |
| `--adapters` | `-a` | Custom file with custom adapter sequences | Optional | Flags non-standard sequencing adapters |

### ⚠️ 1.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `java.lang.OutOfMemoryError: Java heap space` | FastQC exceeded default Java memory (512MB) when analyzing huge files or long reads. | Increase FastQC Java heap memory via environment: `export _JAVA_OPTIONS="-Xmx2048m"` or reduce `-t`. |
| `gzip: unexpected end of file` / `invalid compressed data` | FASTQ file is truncated or corrupted due to incomplete transfer/download. | Verify file integrity using `gzip -t <file.fastq.gz>`. Re-download or re-copy original sequencing files. |
| `Too many open files` | Running FastQC on hundreds of cells simultaneously exceeds system file descriptor limit. | Use `xargs -n 20` to batch files (as implemented in `scripts/qc_data_00.sh`) or increase `ulimit -n 4096`. |
| `Failed to process file ... / permission denied` | Output folder is read-only or does not exist. | Run `mkdir -p <outdir>` and check write permissions with `chmod u+w <outdir>`. |

### 📊 1.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `sample_R1_fastqc.html`: Interactive browser report.
  - `sample_R1_fastqc.zip`: Archive containing raw tabular data (`fastqc_data.txt`).
- **Key Metrics & Quality Benchmarks**:
  - **Per Base Sequence Quality**: Phred scores should be > 28 (green zone). Quality drop-off at the 3' end is typical for Illumina sequencing.
  - **Per Sequence GC Content**: Should display a smooth bell curve matching target organism GC content (~40-50% for human/mouse). Multimodal curves suggest bacterial contamination or primer dimers.
  - **Adapter Content**: Raw libraries often show Nextera/Illumina adapter ramp-up after 50-75bp. This must drop to **0.0%** after fastp trimming.

---

## 2. fastp (Adapter, Poly-G/X, and Quality Trimming)

### 📌 2.1 Purpose
fastp is an ultra-fast, all-in-one FASTQ preprocessor. It performs adapter trimming, poly-G tail removal (NovaSeq artifact), poly-X homopolymer trimming, and quality sliding window filtering.

### 💻 2.2 How to Run
```bash
fastp \
    --in1 data/raw/sample_01_R1.fastq.gz \
    --in2 data/raw/sample_01_R2.fastq.gz \
    --out1 data/clean/sample_01_R1.clean.fastq.gz \
    --out2 data/clean/sample_01_R2.clean.fastq.gz \
    --detect_adapter_for_pe \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x --poly_x_min_len 10 \
    --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
    --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --length_required 35 \
    --thread 8 \
    --json reports/fastp/sample_01_fastp.json \
    --html reports/fastp/sample_01_fastp.html
```

### 🏷️ 2.3 Flag Breakdown
| Flag | Description | Recommended Value | Impact |
|---|---|---|---|
| `--in1`, `--in2` | Input Read 1 and Read 2 FASTQ files | Real paths | Must be paired-end files |
| `--out1`, `--out2` | Output clean Read 1 and Read 2 FASTQ files | `data/clean/*` | Gzip compressed output |
| `--detect_adapter_for_pe` | Auto-detects adapter sequence by read overlap | Enabled | Identifies any Illumina/Nextera adapter automatically |
| `--adapter_sequence` | Explicit Read 1 adapter sequence | `CTGTCTCTTATACACATCT` | Nextera transposase sequence for Smart-seq2 |
| `--adapter_sequence_r2`| Explicit Read 2 adapter sequence | `CTGTCTCTTATACACATCT` | Reverse Nextera sequence |
| `--trim_poly_g` | Detects and trims poly-G tails | Enabled | Essential for NovaSeq/NextSeq 2-color chemistry |
| `--poly_g_min_len` | Minimum length of poly-G tail to trigger trimming | `10` | Prevents false-positive trimming of natural G stretches |
| `--trim_poly_x` | Detects and trims poly-A/T/C homopolymer tails | Enabled | Trims poly-A tails from cDNA synthesis |
| `--poly_x_min_len` | Minimum length of poly-X tail | `10` | Standard threshold |
| `--cut_front` | 5' sliding-window quality trimming | Enabled | Trims low-quality bases at read start |
| `--cut_tail` | 3' sliding-window quality trimming | Enabled | Trims low-quality bases at read ends |
| `--cut_window_size` | Size of the sliding window | `4` | Standard sliding window size |
| `--cut_mean_quality`| Minimum mean Phred quality required in window | `20` (Q20 = 99% accuracy) | Drops window if mean quality < 20 |
| `-q, --qualified_quality_phred` | Phred score required for a base to be qualified | `20` | Standard quality threshold |
| `-u, --unqualified_percent_limit` | Maximum % of unqualified bases before dropping read | `30` | Discards reads with > 30% low-quality bases |
| `-l, --length_required` | Minimum read length to keep after trimming | `35` (or `50` for 150bp runs) | Prevents short ambiguous reads from reaching aligner |
| `--thread` | Number of worker threads | `4 - 8` | Highly efficient multithreading |
| `-j, --json` | JSON output report path | `reports/fastp/*.json` | Parsed directly by MultiQC |
| `-h, --html` | Interactive visual HTML report | `reports/fastp/*.html` | Quality plots before vs after |

### ⚠️ 2.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `Failed to open input file` | Incorrect FASTQ file path or broken symlink. | Ensure path exists and has read permissions. |
| `R1 and R2 have different number of reads` | Pairing mismatch or one mate file was interrupted during download. | Check line counts: `zcat R1.fastq.gz \| wc -l` must equal `zcat R2.fastq.gz \| wc -l`. |
| `> 90% of reads filtered out (too short)` | Read length requirement (`-l`) is set higher than actual sequencer cycle length (e.g. read length was 50bp but `-l` was 75bp). | Inspect raw read length in FastQC and set `-l 35` or `-l 25`. |
| `Excessive poly-G trimming warning` | Sequencing cDNA libraries with long natural G-stretches or 2-color NovaSeq artifact. | Increase `--poly_g_min_len 15` if natural poly-G tracts are expected. |

### 📊 2.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `sample_R1.clean.fastq.gz`, `sample_R2.clean.fastq.gz`: Cleaned paired-end reads ready for alignment.
  - `sample_fastp.json`: Machine-readable trimming statistics.
  - `sample_fastp.html`: Interactive visual dashboard.
- **Key Metrics & Quality Benchmarks**:
  - **Filtering Rate**: > 85-95% of reads should pass filtering.
  - **Q30 After Filtering**: Should reach **> 90-95%**.
  - **Insert Size Peak**: Typically 150-400 bp for Smart-seq2 libraries.

---

## 3. STAR (Splice-Aware Reference Alignment)

### 📌 3.1 Purpose
STAR (Spliced Transcripts Alignment to a Reference) is an ultrafast RNA-seq aligner capable of discovering non-canonical splices and chimeric transcripts using a suffix array indexing strategy.

### 💻 3.2 How to Run
#### Indexing (Run Once):
```bash
STAR \
    --runMode genomeGenerate \
    --genomeDir data/reference/star_index \
    --genomeFastaFiles data/reference/genome.fa \
    --sjdbGTFfile data/reference/genes.gtf \
    --sjdbOverhang 100 \
    --runThreadN 16
```

#### Alignment:
```bash
STAR \
    --runThreadN 16 \
    --genomeDir data/reference/star_index \
    --readFilesIn data/clean/sample_01_R1.clean.fastq.gz data/clean/sample_01_R2.clean.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM NM MD jM jI XS \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --limitBAMsortRAM 31000000000 \
    --quantMode GeneCounts \
    --outFileNamePrefix data/aligned/sample_01_
```

### 🏷️ 3.3 Flag Breakdown
| Flag | Description | Recommended Value | Impact |
|---|---|---|---|
| `--runMode` | Execution mode | `genomeGenerate` or `alignReads` | Default is `alignReads` |
| `--genomeDir` | Path to STAR index directory | Directory path | Must contain pre-built STAR index |
| `--genomeFastaFiles` | Genomic FASTA for index generation | `genome.fa` | Primary assembly FASTA |
| `--sjdbGTFfile` | Gene transfer format (GTF) annotation | `genes.gtf` | Required for splice junction database |
| `--sjdbOverhang` | Splice junction overhang | ReadLength - 1 (e.g. `100`) | Standard is 100 for 100-150bp reads |
| `--readFilesIn` | Input paired-end FASTQ reads | `R1.clean.fastq.gz R2.clean.fastq.gz` | Space-separated paths |
| `--readFilesCommand` | Command to decompress input files | `zcat` (Linux) or `gzcat` (macOS) | Allows reading gzipped FASTQ directly |
| `--outSAMtype` | Output alignment format | `BAM SortedByCoordinate` | Emits coordinate-sorted BAM directly (no samtools sort needed) |
| `--outSAMunmapped` | Handling unmapped reads | `Within` | Keeps unmapped reads in BAM (preserves library total counts) |
| `--outSAMattributes` | SAM tags included in BAM | `NH HI AS nM NM MD jM jI XS` | Essential tags for mapping quality and IGV |
| `--outFilterType` | Splice junction filtering mode | `BySJout` | Filters out spurious junctions |
| `--outFilterMultimapNmax` | Maximum loci for multimapping reads | `20` | Reads mapping to >20 loci are marked unmapped |
| `--outFilterMismatchNmax` | Maximum mismatches allowed per read | `10` | Discards low-similarity mappings |
| `--alignIntronMin` | Minimum intron size | `20` | Canonical mammalian lower intron limit |
| `--alignIntronMax` | Maximum intron size | `1000000` (1 Mb) | Upper intron bound |
| `--alignMatesGapMax` | Maximum genomic distance between mates | `1000000` (1 Mb) | Accommodates large mammalian introns |
| `--limitBAMsortRAM` | RAM buffer for coordinate sorting | `31000000000` (~31 GB) | **CRITICAL**: Prevents OOM crash when sorting human chromosomes |
| `--quantMode` | Internal quantification mode | `GeneCounts` | Emits raw read-to-gene counts table directly from STAR |
| `--outFileNamePrefix` | Output files prefix | `data/aligned/<sample>_` | Standard prefix ending with `_` |

### ⚠️ 3.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `EXITING: FATAL INPUT ERROR: could not open genome file ... Genome` | The `--genomeDir` path is incorrect, or index building was incomplete. | Run `bash scripts/setup_reference.sh` to generate the complete STAR index. |
| `BAMsort: limitBAMsortRAM reached` / `Killed` / Exit code 137 | System ran out of memory while sorting BAM coordinate blocks. | Increase `--limitBAMsortRAM` if RAM allows, or reduce `--runThreadN` to lower concurrent memory usage. |
| `EXITING because of FATAL ERROR: genomeSAindexNbases 6 is too small` | Using toy `--genomeSAindexNbases 6` on a real human/mouse genome. | Remove `--genomeSAindexNbases 6` (STAR defaults to 14 for human/mouse). |
| `Segmentation fault (core dumped)` during index loading | Machine has less than 32GB RAM for human GRCh38 index. | Human STAR index requires ~30-32GB RAM. Ensure machine has $\ge 32$GB RAM or build index with `--genomeSAsparseD 2` to reduce memory to ~16GB. |
| `Uniquely mapped reads % < 40%` | Species mismatch (e.g. mouse reads mapped to human index) or severe DNA contamination. | Verify sequencing source and reference genome matching. |

### 📊 3.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `<sample>_Aligned.sortedByCoord.out.bam`: Coordinate-sorted BAM alignment file.
  - `<sample>_Log.final.out`: Summary statistics of alignment.
  - `<sample>_ReadsPerGene.out.tab`: Raw counts per gene for unstranded and stranded libraries.
  - `<sample>_SJ.out.tab`: High-confidence splice junctions.
- **Key Metrics & Quality Benchmarks**:
  - **Uniquely Mapped Reads %**: **> 70-85%** (Excellent). < 50% indicates degradation or wrong reference.
  - **% of Reads Mapped to Multiple Loci**: Typically 5-15% (normal for repetitive/homologous genes).
  - **Mismatch Rate per Base**: Should be **< 0.5% - 1.0%**.

---

## 4. samtools (BAM Sorting & Coordinate Indexing)

### 📌 4.1 Purpose
samtools provides utilities for manipulating alignments in the SAM/BAM format, including indexing, coordinate sorting, flag extraction, and coverage statistics.

### 💻 4.2 How to Run
```bash
# Index a coordinate-sorted BAM:
samtools index -@ 8 data/aligned/sample_01_Aligned.sortedByCoord.out.bam

# Alignment statistics check:
samtools flagstat data/aligned/sample_01_Aligned.sortedByCoord.out.bam

# Count mapped reads:
samtools view -c -F 4 data/aligned/sample_01_Aligned.sortedByCoord.out.bam
```

### 🏷️ 4.3 Flag Breakdown
| Command / Flag | Description | Recommended Value | Impact |
|---|---|---|---|
| `samtools index` | Creates `.bai` coordinate index | BAM path | Allows instantaneous coordinate querying |
| `-@` | Number of indexing threads | `4 - 8` | Accelerates index generation |
| `samtools flagstat` | Prints overall alignment summary | BAM path | Reports total, mapped, paired, duplicated reads |
| `-F 4` | Exclude unmapped reads | Flag filter | Counts only mapped reads |

### ⚠️ 4.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `[E::hts_idx_push] NO_COOR reads not in a single block at the end` | The BAM file is queryname-sorted or unsorted, not coordinate-sorted. | Run `samtools sort -o sorted.bam input.bam` before indexing. |
| `[E::hts_open_format] Failed to open file ...: No such file or directory` | BAM path is wrong or previous STAR alignment failed. | Verify that STAR alignment completed without errors. |
| `[W::bam_hdr_read] EOF marker is absent` | BAM file was truncated (e.g. process was killed midway). | Alignment was incomplete. Re-run alignment stage. |

### 📊 4.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `<sample>_Aligned.sortedByCoord.out.bam.bai`: Binary index file allowing fast random access.
- **Key Metrics**:
  - `samtools flagstat`: Should show `100% paired` and `> 80% mapped`.

---

## 5. Subread featureCounts (Gene Quantification)

### 📌 5.1 Purpose
featureCounts is an efficient read summarization program that quantifies the number of reads or read pairs overlapping genomic features (such as exons and genes) defined in a GTF annotation file.

### 💻 5.2 How to Run
```bash
featureCounts \
    -T 16 \
    -p \
    --countReadPairs \
    -t exon \
    -g gene_id \
    -a data/reference/genes.gtf \
    -o data/counts/featurecounts_raw.txt \
    -s 0 \
    -Q 10 \
    data/aligned/*.sortedByCoord.out.bam
```

### 🏷️ 5.3 Flag Breakdown
| Flag | Description | Recommended Value | Impact |
|---|---|---|---|
| `-T` | Number of parallel threads | `8 - 16` | Fast multi-sample quantification |
| `-p` | Paired-end mode | Enabled | Treats mate pairs as fragments |
| `--countReadPairs` | Count fragments instead of individual reads | Enabled | **CRITICAL**: Prevents counting both ends of a pair twice |
| `-t` | Feature type in GTF (column 3) | `exon` | Quantifies reads falling on annotated exons |
| `-g` | Attribute in GTF (column 9) to group by | `gene_id` | Aggregates all exons belonging to the same gene |
| `-a` | Path to GTF annotation file | `data/reference/genes.gtf` | Must match chromosome names of BAM files |
| `-o` | Output raw counts text file | `data/counts/featurecounts_raw.txt` | Tab-delimited quantification table |
| `-s` | Strand specificity | `0` (Unstranded for Smart-seq2) | `0` = unstranded, `1` = stranded, `2` = reverse stranded |
| `-Q` | Minimum mapping quality (MAPQ) threshold | `10` | Filters out low-confidence/multimapping reads |
| `--primary` | Only count primary alignments | Recommended | Avoids counting secondary alignments |
| `-B` | Require both ends of pair to be aligned | Optional | Stricter paired-end requirement |
| `-C` | Exclude chimeric fragments | Optional | Ignores pairs where mates align to different chromosomes |

### ⚠️ 5.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `Assigned: 0 (0.0%)` across all samples! | **Chromosome naming mismatch** between BAM and GTF! e.g. BAM has `chr1, chr2` while GTF has `1, 2`. | Check BAM header: `samtools view -H sample.bam \| grep "@SQ"`. Check GTF: `grep -v "^#" genes.gtf \| head -n 5`. Align naming or use `setup_reference.sh`. |
| `Assigned: 0%` with non-zero alignment | Incorrect strand specificity (`-s 1` or `-s 2` used on unstranded library). | For Smart-seq2, always use `-s 0` (unstranded). |
| `Failed to open annotation file ...` | GTF file does not exist or is gzipped (`.gtf.gz`). | Decompress GTF: `gzip -d genes.gtf.gz`. featureCounts requires uncompressed GTF. |
| `Unknown option --countReadPairs` | Subread version is outdated (< 1.6). | Update subread via pixi (`pixi install subread>=2.0.6`). |

### 📊 5.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `featurecounts_raw.txt`: Full table with columns `Geneid, Chr, Start, End, Strand, Length, <sample_01>...`.
  - `featurecounts_raw.txt.summary`: Summary breakdown table.
  - `gene_cell_count_matrix.tsv`: Clean Gene $\times$ Cell expression matrix formatted via AWK (Geneid in column 1, cell names across header, integer counts).
- **Key Metrics & Quality Benchmarks**:
  - **Assigned Rate**: **> 60-80%** of alignments should be `Assigned`.
  - **Unassigned_NoFeatures**: Normal range 10-25%. > 40% indicates high genomic DNA contamination or intron retention.
  - **Unassigned_Ambiguity**: Typically < 5-10% (overlapping gene annotations).

---

## 6. MultiQC (Comprehensive Metric Aggregation)

### 📌 6.1 Purpose
MultiQC searches a given directory for analysis logs across multiple bioinformatics tools and compiles them into a single, cohesive, interactive HTML report.

### 💻 6.2 How to Run
```bash
multiqc \
    --outdir reports/multiqc \
    --filename single_cell_multiqc_report.html \
    --force \
    --interactive \
    reports/qc_raw reports/fastp data/aligned data/counts
```

### 🏷️ 6.3 Flag Breakdown
| Flag | Description | Recommended Value | Impact |
|---|---|---|---|
| `--outdir` (`-o`) | Directory to write the report | `reports/multiqc` | Output folder |
| `--filename` (`-n`) | Custom name for the HTML report file | `single_cell_multiqc_report.html` | Custom report name |
| `--force` (`-f`) | Overwrites existing report files | Enabled | Prevents report naming collisions |
| `--interactive` | Forces all plots to be interactive | Enabled | Allows zooming and filtering directly in report |
| `--dirs` | Prepends directory names to sample names | Optional | Prevents collisions if duplicate sample names exist |
| `--config` | Custom configuration YAML | Optional | Custom colors, logos, or section ordering |
| `<dirs...>` | Input search directories | Space-separated list | MultiQC parses all logs recursively |

### ⚠️ 6.4 Common Errors & Troubleshooting ("เจอ error อะไร")
| Error Message / Symptom | Root Cause | Concrete Solution |
|---|---|---|
| `No analysis results found` | Input directories do not contain recognizable log files, or path is incorrect. | Pass the exact folders containing logs: `reports/qc_raw`, `reports/fastp`, `data/aligned`, `data/counts`. |
| `MemoryError / Browser freezes` | Visualizing thousands of single cells at once in interactive mode. | Remove `--interactive` for large cohorts (> 1,000 cells) to enable flat rendering. |
| `Duplicate sample names warning` | Sample names conflict across different subdirectories. | Add `--dirs` flag to qualify sample names with their directory path. |

### 📊 6.5 Expected Outputs & Metric Interpretation
- **Output Files**:
  - `single_cell_multiqc_report.html`: Interactive self-contained HTML report.
  - `single_cell_multiqc_report_data/`: Directory containing parsed stats (`multiqc_data.json`, `multiqc_general_stats.txt`, `multiqc_star.txt`, `multiqc_fastp.txt`).
- **Key Modules Displayed**:
  1. **General Statistics Table**: Summary of raw reads, trimming loss, mapping rates, and feature assignment per cell.
  2. **FastQC**: Per-base quality scores and sequence duplication.
  3. **fastp**: Filtering results and insert size distributions.
  4. **STAR**: Uniquely mapped reads vs multimapping breakdown.
  5. **featureCounts**: Exonic assigned read percentages per cell.
