#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/quant_03.sh
# Purpose: Gene-Level Quantification via Subread featureCounts (Production Ready)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${ALIGNED_DIR}}"
OUT_DIR="${2:-${COUNTS_DIR}}"
ANNOTATION_GTF="${3:-${GTF_FILE}}"
CORE_THREADS="${4:-${THREADS}}"

log_step "Starting Stage 03: Production Subread featureCounts Quantification"
log_info "Aligned BAM directory: ${IN_DIR}"
log_info "Counts directory:      ${OUT_DIR}"
log_info "GTF Annotation:        ${ANNOTATION_GTF}"
log_info "Threads:               ${CORE_THREADS}"
log_info "Feature / Attribute:   ${FEATURE_TYPE} / ${ATTRIBUTE_TYPE}"
log_info "Strandedness:          ${STRANDEDNESS}"
log_info "Min MAPQ:              ${MIN_MAPQ}"

# 1. Validate GTF file
if [ ! -f "${ANNOTATION_GTF}" ]; then
    log_error "GTF annotation file not found: ${ANNOTATION_GTF}!"
    log_info "Please download/provide a valid GTF file (e.g. via scripts/setup_reference.sh)"
    exit 1
fi

mkdir -p "${OUT_DIR}"

# 2. Collect all BAM files
mapfile -t BAM_FILES < <(find "${IN_DIR}" -maxdepth 2 -type f -name "*_Aligned.sortedByCoord.out.bam" | sort)

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    log_error "No sorted BAM files found in ${IN_DIR}! Please run scripts/align_02.sh first."
    exit 1
fi

log_info "Quantifying expression across ${#BAM_FILES[@]} BAM libraries..."

RAW_COUNTS="${OUT_DIR}/featurecounts_raw.txt"
CLEAN_MATRIX="${OUT_DIR}/gene_cell_count_matrix.tsv"

# Run featureCounts across all BAM files simultaneously
featureCounts \
    -T "${CORE_THREADS}" \
    -p \
    --countReadPairs \
    -t "${FEATURE_TYPE}" \
    -g "${ATTRIBUTE_TYPE}" \
    -a "${ANNOTATION_GTF}" \
    -o "${RAW_COUNTS}" \
    -s "${STRANDEDNESS}" \
    -Q "${MIN_MAPQ}" \
    "${BAM_FILES[@]}"

log_ok "Raw featureCounts summary and table written to ${RAW_COUNTS}."

# 3. Format clean Gene x Cell Expression Matrix using pure AWK
log_info "Formatting clean Gene x Cell expression matrix..."

awk '
BEGIN { FS="\t"; OFS="\t" }
/^#/ { next }
NR==2 {
    # Header row: print Geneid followed by simplified sample names
    printf "%s", $1
    for (i=7; i<=NF; i++) {
        n = split($i, parts, "/")
        fname = parts[n]
        sub(/_Aligned\.sortedByCoord\.out\.bam$/, "", fname)
        printf "\t%s", fname
    }
    printf "\n"
    next
}
{
    # Data rows: Geneid and count values
    printf "%s", $1
    for (i=7; i<=NF; i++) {
        printf "\t%s", $i
    }
    printf "\n"
}
' "${RAW_COUNTS}" > "${CLEAN_MATRIX}"

log_ok "Gene-cell expression matrix created successfully at: ${CLEAN_MATRIX}"
