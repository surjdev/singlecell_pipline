#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/align_02.sh
# Purpose: Splice-Aware Genome Alignment via STAR & BAM Indexing (Production Ready)
# ==============================================================================
# 🏷️ Flags & Usage:
#   --runThreadN <N>             : CPU threads for parallel alignment and sorting
#   --genomeDir <dir>            : Path to pre-built STAR suffix array index
#   --readFilesIn <R1> <R2>      : Space-separated paired-end cleaned FASTQ files
#   --readFilesCommand zcat      : On-the-fly decompression command for gzipped FASTQs
#   --outSAMtype BAM SortedByCoordinate : Direct coordinate-sorted BAM output
#   --outSAMunmapped Within      : Retain unmapped reads in BAM (preserves library total counts)
#   --outSAMattributes NH HI AS nM NM MD jM jI XS : Essential tags for MAPQ, IGV, and single-cell
#   --outFilterType BySJout      : Filter spurious splice junctions
#   --outFilterMultimapNmax 20   : Discard reads mapping to >20 loci
#   --outFilterMismatchNmax 10   : Max mismatch allowance per read
#   --alignIntronMin 20          : Lower bound for mammalian splice junctions
#   --alignIntronMax 1000000     : Upper bound for mammalian splice junctions (1 Mb)
#   --alignMatesGapMax 1000000   : Max distance between mate pairs across introns (1 Mb)
#   --limitBAMsortRAM <bytes>    : Sorting buffer memory (~31GB). Prevents OOM crash
#   --quantMode GeneCounts       : Quantify raw read-to-gene counts directly inside STAR
#   --outFileNamePrefix <prefix> : Output files prefix ending with '_'
#   samtools index -@ <threads>  : Create BAM coordinate index (.bai)
#
# ⚠️ Common Errors & Fixes ("เจอ error อะไร"):
#   1. "BAMsort: limitBAMsortRAM reached / Killed / Exit code 137":
#      - Cause: Insufficient sorting buffer or system RAM during chromosome sorting.
#      - Fix: Increase --limitBAMsortRAM or reduce --runThreadN to lower concurrent RAM.
#   2. "FATAL INPUT ERROR: could not open genome file ... Genome":
#      - Cause: Incomplete or missing STAR index.
#      - Fix: Run scripts/setup_reference.sh to download reference and build full index.
#   3. "Segmentation fault (core dumped)":
#      - Cause: System has < 32GB RAM for human GRCh38 genome index.
#      - Fix: Ensure machine has >= 32GB RAM or use --genomeSAsparseD 2.
#   4. "Uniquely mapped reads % < 40%":
#      - Cause: Reference species mismatch (e.g. mouse reads aligned to human index).
#      - Fix: Check sample species and verify reference index genome.
#
# 📊 Expected Outputs ("ผลลัพธ์เป็นอย่างไร"):
#   - data/aligned/*_Aligned.sortedByCoord.out.bam : Coordinate-sorted BAM alignment file
#   - data/aligned/*_Aligned.sortedByCoord.out.bam.bai : samtools coordinate index
#   - data/aligned/*_Log.final.out               : Alignment rate and error summary
#   - data/aligned/*_ReadsPerGene.out.tab        : STAR gene counts table
#   - data/aligned/*_SJ.out.tab                  : High-confidence splice junctions
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${CLEAN_DIR}}"
OUT_DIR="${2:-${ALIGNED_DIR}}"
INDEX_DIR="${3:-${STAR_INDEX}}"
CORE_THREADS="${4:-${THREADS}}"

[[ "$CORE_THREADS" =~ ^[1-9][0-9]*$ ]] || { log_error "threads must be positive"; exit 2; }
[[ -d "$IN_DIR" ]] || { log_error "Input directory does not exist: $IN_DIR"; exit 2; }

log_step "Starting Stage 02: Production STAR Splice-Aware Alignment"
log_info "Clean FASTQ directory: ${IN_DIR}"
log_info "Aligned BAM directory: ${OUT_DIR}"
log_info "STAR Reference Index:  ${INDEX_DIR}"
log_info "Threads:               ${CORE_THREADS}"
log_info "RAM limit for sorting: $((STAR_RAM_LIMIT / 1024 / 1024 / 1024)) GB"

# 1. Validate STAR Reference Index
if [ ! -f "${INDEX_DIR}/SA" ] || [ ! -f "${INDEX_DIR}/Genome" ]; then
    log_error "STAR index not found or incomplete in ${INDEX_DIR}!"
    log_info "Please generate or download a production STAR index first, e.g.:"
    log_info "  bash scripts/setup_reference.sh --species human --threads ${CORE_THREADS}"
    log_info "Or point STAR_INDEX to an existing index in config.env or environment."
    exit 1
fi

mkdir -p "${OUT_DIR}"

# 2. Identify all cleaned paired-end Read 1 files
mapfile -t CLEAN_R1 < <(find -L "${IN_DIR}" -maxdepth 2 -type f \( \
    -name "*_R1.clean.fastq.gz" -o \
    -name "*_1.clean.fastq.gz"  -o \
    -name "*_R1.clean.fq.gz"    -o \
    -name "*_1.clean.fq.gz"        \
\) | sort)

if [ ${#CLEAN_R1[@]} -eq 0 ]; then
    log_error "No cleaned paired-end FASTQ files found in ${IN_DIR}!"
    log_info "Please run scripts/fastp_01.sh first."
    exit 1
fi

log_info "Aligning ${#CLEAN_R1[@]} clean single-cell libraries..."

# 3. Align each library
ALIGNED_COUNT=0
for R1 in "${CLEAN_R1[@]}"; do
    DIRNAME=$(dirname "${R1}")
    FILENAME=$(basename "${R1}")

    if [[ "${FILENAME}" =~ ^(.*)_R1\.clean\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_R2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_1\.clean\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_R1\.clean\.fq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_R2.clean.fq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_1\.clean\.fq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_2.clean.fq.gz"
    else
        log_warn "Unrecognized naming format for ${FILENAME}, skipping."
        continue
    fi

    if [ ! -f "${R2}" ]; then
        log_error "Missing paired-end R2 mate for ${R1} (expected: ${R2})!"
        exit 1
    fi

    PREFIX="${OUT_DIR}/${SAMPLE_ID}_"
    SORTED_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"

    log_info "Aligning library: ${SAMPLE_ID}"

    STAR \
        --runThreadN "${CORE_THREADS}" \
        --genomeDir "${INDEX_DIR}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD jM jI XS \
        --quantMode GeneCounts \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --limitBAMsortRAM "${STAR_RAM_LIMIT}" \
        --outFileNamePrefix "${PREFIX}"

    # Index sorted BAM with samtools
    log_info "Indexing BAM: ${SORTED_BAM}"
    samtools index -@ "${CORE_THREADS}" "${SORTED_BAM}"

    ALIGNED_COUNT=$((ALIGNED_COUNT + 1))
done

log_ok "Completed STAR alignment and BAM indexing for ${ALIGNED_COUNT} libraries."
log_info "BAM alignments saved to: ${OUT_DIR}"
