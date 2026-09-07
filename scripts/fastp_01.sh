#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/fastp_01.sh
# Purpose: Adapter, Poly-G/X, and Quality Trimming via fastp (Production Ready)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${RAW_DIR}}"
OUT_DIR="${2:-${CLEAN_DIR}}"
REPORT_DIR="${3:-${FASTP_DIR}}"
CORE_THREADS="${4:-${THREADS}}"

log_step "Starting Stage 01: Production fastp Quality & Adapter Trimming"
log_info "Input directory:       ${IN_DIR}"
log_info "Clean FASTQ directory: ${OUT_DIR}"
log_info "Report directory:      ${REPORT_DIR}"
log_info "Threads:               ${CORE_THREADS}"
log_info "Min read length:       ${MIN_READ_LENGTH} bp"
log_info "Min Phred quality:     Q${MIN_QUALITY}"

mkdir -p "${OUT_DIR}" "${REPORT_DIR}"

# Find all Read 1 files with diverse real-world Illumina naming patterns
mapfile -t R1_FILES < <(find "${IN_DIR}" -maxdepth 2 -type f \( \
    -name "*_R1_001.fastq.gz" -o \
    -name "*_R1.fastq.gz"     -o \
    -name "*.R1.fastq.gz"     -o \
    -name "*_1.fastq.gz"      -o \
    -name "*_1.fq.gz"         -o \
    -name "*_R1.fq.gz"           \
\) | sort)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    log_error "No paired-end R1 FASTQ files found in ${IN_DIR}!"
    log_info "Supported filename patterns: *_R1.fastq.gz, *_R1_001.fastq.gz, *_1.fastq.gz, *_1.fq.gz"
    exit 1
fi

log_info "Identified ${#R1_FILES[@]} paired-end libraries to process."

# Process each library pair
PROCESSED_COUNT=0
for R1 in "${R1_FILES[@]}"; do
    DIRNAME=$(dirname "${R1}")
    FILENAME=$(basename "${R1}")
    
    # Robust sample name extraction and mate matching
    if [[ "${FILENAME}" =~ ^(.*)_R1_001\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_R2_001.fastq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_R1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_R2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_R1\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_R2.fastq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_R1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_R2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)\.R1\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}.R2.fastq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_R1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_R2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_1\.fastq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_2.fastq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_1\.fq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_2.fq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_2.clean.fastq.gz"
    elif [[ "${FILENAME}" =~ ^(.*)_R1\.fq\.gz$ ]]; then
        SAMPLE_ID="${BASH_REMATCH[1]}"
        R2="${DIRNAME}/${SAMPLE_ID}_R2.fq.gz"
        OUT_R1="${OUT_DIR}/${SAMPLE_ID}_R1.clean.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE_ID}_R2.clean.fastq.gz"
    else
        log_warn "Could not deduce pair pattern for: ${FILENAME}, skipping."
        continue
    fi

    if [ ! -f "${R2}" ]; then
        log_error "Missing paired-end mate for ${R1} (expected: ${R2})!"
        exit 1
    fi

    JSON_REPORT="${REPORT_DIR}/${SAMPLE_ID}_fastp.json"
    HTML_REPORT="${REPORT_DIR}/${SAMPLE_ID}_fastp.html"

    log_info "Processing sample: ${SAMPLE_ID}"

    # Build adapter arguments
    ADAPTER_ARGS=("--detect_adapter_for_pe")
    if [ "${ADAPTER_FWD}" != "auto" ] && [ -n "${ADAPTER_FWD}" ]; then
        ADAPTER_ARGS+=("--adapter_sequence" "${ADAPTER_FWD}")
        if [ "${ADAPTER_REV}" != "auto" ] && [ -n "${ADAPTER_REV}" ]; then
            ADAPTER_ARGS+=("--adapter_sequence_r2" "${ADAPTER_REV}")
        fi
    fi

    fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${OUT_R1}" \
        --out2 "${OUT_R2}" \
        "${ADAPTER_ARGS[@]}" \
        --trim_poly_g \
        --poly_g_min_len 10 \
        --trim_poly_x \
        --poly_x_min_len 10 \
        --cut_front \
        --cut_front_window_size 4 \
        --cut_front_mean_quality "${MIN_QUALITY}" \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality "${MIN_QUALITY}" \
        --qualified_quality_phred "${MIN_QUALITY}" \
        --unqualified_percent_limit "${MAX_UNQUALIFIED_PCT}" \
        --length_required "${MIN_READ_LENGTH}" \
        --thread "${CORE_THREADS}" \
        --json "${JSON_REPORT}" \
        --html "${HTML_REPORT}"

    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
done

log_ok "Completed fastp trimming for ${PROCESSED_COUNT} libraries."
log_info "Clean FASTQ outputs saved to: ${OUT_DIR}"
log_info "Quality reports saved to:     ${REPORT_DIR}"
