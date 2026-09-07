#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/qc_data_00.sh
# Purpose: Raw FASTQ Read Quality Control via FastQC (Production Ready)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${RAW_DIR}}"
OUT_DIR="${2:-${QC_RAW_DIR}}"
CORE_THREADS="${3:-${THREADS}}"

log_step "Starting Stage 00: Raw Read Quality Control (FastQC)"
log_info "Input directory:  ${IN_DIR}"
log_info "Output directory: ${OUT_DIR}"
log_info "Threads:          ${CORE_THREADS}"

mkdir -p "${OUT_DIR}"

# Collect all FASTQ variations (.fastq.gz, .fq.gz, .fastq, .fq)
mapfile -t FASTQ_FILES < <(find "${IN_DIR}" -maxdepth 2 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) | sort)

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    log_error "No FASTQ files found in ${IN_DIR}!"
    log_info "Please place paired-end FASTQ files in ${IN_DIR} (e.g. sample_01_R1.fastq.gz, sample_01_R2.fastq.gz)"
    exit 1
fi

log_info "Found ${#FASTQ_FILES[@]} FASTQ files to analyze."

# Execute FastQC with parallel threads and batching
printf '%s\n' "${FASTQ_FILES[@]}" | xargs -n 20 -P 1 fastqc \
    --outdir "${OUT_DIR}" \
    --threads "${CORE_THREADS}" \
    --noextract \
    --quiet

log_ok "FastQC completed successfully. Reports saved to ${OUT_DIR}"
