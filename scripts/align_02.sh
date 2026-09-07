#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/align_02.sh
# Purpose: Splice-Aware Genome Alignment via STAR & BAM Indexing (Production Ready)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${CLEAN_DIR}}"
OUT_DIR="${2:-${ALIGNED_DIR}}"
INDEX_DIR="${3:-${STAR_INDEX}}"
CORE_THREADS="${4:-${THREADS}}"

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
mapfile -t CLEAN_R1 < <(find "${IN_DIR}" -maxdepth 2 -type f \( \
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
