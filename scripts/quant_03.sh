#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/quant_03.sh
# Purpose: Gene-Level Quantification via Subread featureCounts (Production Ready)
# ==============================================================================
# 🏷️ Flags & Usage:
#   -T <threads>          : CPU threads for parallel BAM parsing
#   -p                    : Paired-end mode (treats mate pairs as single fragments)
#   --countReadPairs      : Count fragments rather than individual reads twice
#   -t exon               : Target feature type in column 3 of GTF (default: exon)
#   -g gene_id            : Meta-feature grouping attribute in column 9 of GTF (default: gene_id)
#   -a <gtf>              : Reference annotation GTF path
#   -o <file>             : Output raw counts file path
#   -s <0|1|2>            : Strand specificity (0 = unstranded Smart-seq2, 1 = sense, 2 = antisense)
#   -Q 10                 : Minimum MAPQ mapping quality threshold
#
# ⚠️ Common Errors & Fixes ("เจอ error อะไร"):
#   1. "Assigned: 0 (0.0%) across all samples!":
#      - Cause: Chromosome naming mismatch between BAM (e.g. chr1, chr2) and GTF (e.g. 1, 2)
#               or inverted strand specificity (-s 1 instead of -s 0).
#      - Fix: Compare BAM header (samtools view -H <bam> | grep "@SQ") with GTF header
#             (grep -v "^#" <gtf> | head -n 5) and align naming prefixes.
#   2. "Failed to open annotation file ...":
#      - Cause: GTF file is gzipped (.gtf.gz). featureCounts requires uncompressed GTF.
#      - Fix: Decompress GTF using: gzip -d genes.gtf.gz
#   3. "Unknown option --countReadPairs":
#      - Cause: Subread version is outdated (< 1.6).
#      - Fix: Update subread via pixi: pixi install subread>=2.0.6
#
# 📊 Expected Outputs ("ผลลัพธ์เป็นอย่างไร"):
#   - data/counts/featurecounts_raw.txt         : Full raw quantification table
#   - data/counts/featurecounts_raw.txt.summary : Read assignment summary (Assigned, Unassigned)
#   - data/counts/gene_cell_count_matrix.tsv    : Clean Gene x Cell matrix formatted via AWK
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

IN_DIR="${1:-${ALIGNED_DIR}}"
OUT_DIR="${2:-${COUNTS_DIR}}"
ANNOTATION_GTF="${3:-${GTF_FILE}}"
CORE_THREADS="${4:-${THREADS}}"

[[ "$CORE_THREADS" =~ ^[1-9][0-9]*$ ]] || { log_error "threads must be positive"; exit 2; }
[[ -d "$IN_DIR" ]] || { log_error "Input directory does not exist: $IN_DIR"; exit 2; }

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
mapfile -t BAM_FILES < <(find -L "${IN_DIR}" -maxdepth 2 -type f -name "*_Aligned.sortedByCoord.out.bam" | sort)

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
$1 == "Geneid" {
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
