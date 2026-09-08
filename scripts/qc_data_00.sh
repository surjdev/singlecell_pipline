#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/qc_data_00.sh
# Purpose: Raw FASTQ Read Quality Control via FastQC (Production Ready)
# ==============================================================================
# 🏷️ Flags & Usage:
#   --outdir <dir>  : Destination directory for HTML and ZIP quality reports
#   --threads <num> : Number of parallel files to process (~250MB RAM/thread)
#   --noextract     : Do not decompress output ZIP archive (saves storage/inodes)
#   --quiet         : Suppress verbose progress messages for batch execution
#
# ⚠️ Common Errors & Fixes ("เจอ error อะไร"):
#   1. "java.lang.OutOfMemoryError: Java heap space":
#      - Cause: Default FastQC Java memory (512MB) exceeded by large files.
#      - Fix: Set: export _JAVA_OPTIONS="-Xmx2048m"
#   2. "gzip: unexpected end of file / invalid compressed data":
#      - Cause: FASTQ file truncated during transfer or incomplete download.
#      - Fix: Test with gzip -t <file> and re-transfer corrupted files.
#   3. "Too many open files":
#      - Cause: FastQC invoked on too many files simultaneously.
#      - Fix: Handled automatically in this script via xargs -0 -n 20 batching.
#
# 📊 Expected Outputs ("ผลลัพธ์เป็นอย่างไร"):
#   - reports/qc_raw/*_fastqc.html : Interactive visual QC dashboard per library
#   - reports/qc_raw/*_fastqc.zip  : Raw data tables (fastqc_data.txt) for MultiQC
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${RAW_DIR}}"
OUT_DIR="${2:-${QC_RAW_DIR}}"
CORE_THREADS="${3:-${THREADS}}"

[[ "$CORE_THREADS" =~ ^[1-9][0-9]*$ ]] || { log_error "threads must be positive"; exit 2; }
[[ -d "$IN_DIR" ]] || { log_error "Input directory does not exist: $IN_DIR"; exit 2; }

log_step "Starting Stage 00: Raw Read Quality Control (FastQC)"
log_info "Input directory:  ${IN_DIR}"
log_info "Output directory: ${OUT_DIR}"
log_info "Threads:          ${CORE_THREADS}"

mkdir -p "${OUT_DIR}"

# Collect all FASTQ variations (.fastq.gz, .fq.gz, .fastq, .fq)
mapfile -t FASTQ_FILES < <(find -L "${IN_DIR}" -maxdepth 2 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) | sort)

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    log_error "No FASTQ files found in ${IN_DIR}!"
    log_info "Please place paired-end FASTQ files in ${IN_DIR} (e.g. sample_01_R1.fastq.gz, sample_01_R2.fastq.gz)"
    exit 1
fi

log_info "Found ${#FASTQ_FILES[@]} FASTQ files to analyze."

# Execute FastQC with parallel threads and batching
printf '%s\0' "${FASTQ_FILES[@]}" | xargs -0 -n 20 -P 1 fastqc \
    --outdir "${OUT_DIR}" \
    --threads "${CORE_THREADS}" \
    --noextract \
    --quiet

log_ok "FastQC completed successfully. Reports saved to ${OUT_DIR}"
