#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/multiqc_04.sh
# Purpose: Comprehensive Multi-Tool QC Reporting via MultiQC (Production Ready)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

OUT_DIR="${1:-${MULTIQC_DIR}}"
shift || true
SEARCH_DIRS=("$@")

if [ ${#SEARCH_DIRS[@]} -eq 0 ]; then
    SEARCH_DIRS=("${QC_RAW_DIR}" "${FASTP_DIR}" "${ALIGNED_DIR}" "${COUNTS_DIR}")
fi

log_step "Starting Stage 04: Production MultiQC Report Aggregation"
log_info "Output directory: ${OUT_DIR}"
log_info "Searching in:     ${SEARCH_DIRS[*]}"

mkdir -p "${OUT_DIR}"

multiqc \
    --outdir "${OUT_DIR}" \
    --filename "single_cell_multiqc_report.html" \
    --force \
    --interactive \
    "${SEARCH_DIRS[@]}"

log_ok "MultiQC report compiled successfully!"
log_info "Open in browser: ${OUT_DIR}/single_cell_multiqc_report.html"
