#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/multiqc_04.sh
# Purpose: Comprehensive Multi-Tool QC Reporting via MultiQC (Production Ready)
# ==============================================================================
# 🏷️ Flags & Usage:
#   --outdir <dir>       : Output directory for report HTML and parsed data tables
#   --filename <name>    : Custom HTML report filename (single_cell_multiqc_report.html)
#   --force              : Overwrite previously generated reports
#   --interactive        : Force interactive dynamic charts for zooming and filtering
#   <dir1> <dir2> ...    : Directories searched recursively for analysis logs
#
# ⚠️ Common Errors & Fixes ("เจอ error อะไร"):
#   1. "No analysis results found":
#      - Cause: Search directories do not contain recognizable log files.
#      - Fix: Explicitly pass directories: reports/qc_raw reports/fastp data/aligned data/counts
#   2. "Duplicate sample names warning":
#      - Cause: Identical filenames across different search subfolders.
#      - Fix: MultiQC handles this automatically or pass --dirs to include directory prefixes.
#
# 📊 Expected Outputs ("ผลลัพธ์เป็นอย่างไร"):
#   - reports/multiqc/single_cell_multiqc_report.html : Unified browser-ready dashboard
#   - reports/multiqc/single_cell_multiqc_report_data/ : Raw tabular data tables (JSON, TSV)
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
    --fullnames \
    --interactive \
    "${SEARCH_DIRS[@]}"

log_ok "MultiQC report compiled successfully!"
log_info "Open in browser: ${OUT_DIR}/single_cell_multiqc_report.html"
