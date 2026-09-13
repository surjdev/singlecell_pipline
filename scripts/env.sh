#!/usr/bin/env bash
# ==============================================================================
# Production Environment & Configuration for Upstream Single-Cell Pipeline
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Optional user configuration file override
CONFIG_FILE="${CONFIG_FILE:-${PROJECT_ROOT}/config.env}"
if [ -f "${CONFIG_FILE}" ]; then
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
fi

# ------------------------------------------------------------------------------
# 1. Directory Structure
# ------------------------------------------------------------------------------
RAW_DIR="${RAW_DIR:-${PROJECT_ROOT}/data/raw}"
CLEAN_DIR="${CLEAN_DIR:-${PROJECT_ROOT}/data/clean}"
ALIGNED_DIR="${ALIGNED_DIR:-${PROJECT_ROOT}/data/aligned}"
COUNTS_DIR="${COUNTS_DIR:-${PROJECT_ROOT}/data/counts}"
REPORTS_DIR="${REPORTS_DIR:-${PROJECT_ROOT}/reports}"
REF_DIR="${REF_DIR:-${PROJECT_ROOT}/data/reference}"

QC_RAW_DIR="${QC_RAW_DIR:-${REPORTS_DIR}/qc_raw}"
QC_CLEAN_DIR="${QC_CLEAN_DIR:-${REPORTS_DIR}/qc_clean}"
FASTP_DIR="${FASTP_DIR:-${REPORTS_DIR}/fastp}"
MULTIQC_DIR="${MULTIQC_DIR:-${REPORTS_DIR}/multiqc}"

# ------------------------------------------------------------------------------
# 2. Reference Genome Paths (Human GRCh38 / Mouse GRCm39 / Custom)
# ------------------------------------------------------------------------------
GENOME_FASTA="${GENOME_FASTA:-${REF_DIR}/genome.fa}"
GTF_FILE="${GTF_FILE:-${REF_DIR}/genes.gtf}"
STAR_INDEX="${STAR_INDEX:-${REF_DIR}/star_index}"

# ------------------------------------------------------------------------------
# 3. Hardware Resources
# ------------------------------------------------------------------------------
# Auto-detect available cores if not specified
if [ -z "${THREADS:-}" ]; then
    if command -v nproc >/dev/null 2>&1; then
        THREADS="$(nproc)"
    else
        THREADS=4
    fi
fi
# Cap default threads to avoid system lockup if not explicitly set
THREADS="${THREADS:-4}"

# Maximum RAM for STAR coordinate sorting (in bytes, default ~30GB for mammalian genomes)
STAR_RAM_LIMIT="${STAR_RAM_LIMIT:-31000000000}"

# ------------------------------------------------------------------------------
# 4. Trimming & Quality Thresholds (fastp)
# ------------------------------------------------------------------------------
MIN_READ_LENGTH="${MIN_READ_LENGTH:-25}"       # Minimum read length (25 bp for 43 bp HiSeq 2000 reads)
MIN_QUALITY="${MIN_QUALITY:-20}"               # Phred score >= Q20 (99% base accuracy)
MAX_UNQUALIFIED_PCT="${MAX_UNQUALIFIED_PCT:-30}" # Max % of bases allowed below Q20
ADAPTER_FWD="${ADAPTER_FWD:-CTGTCTCTTATACACATCT}" # Nextera adapter sequence for Smart-seq2
ADAPTER_REV="${ADAPTER_REV:-CTGTCTCTTATACACATCT}" # Nextera adapter sequence

# ------------------------------------------------------------------------------
# 5. Quantification Options (featureCounts)
# ------------------------------------------------------------------------------
STRANDEDNESS="${STRANDEDNESS:-0}"              # 0 = unstranded (Smart-seq2), 1 = stranded, 2 = reverse
FEATURE_TYPE="${FEATURE_TYPE:-exon}"           # GTF feature (column 3)
ATTRIBUTE_TYPE="${ATTRIBUTE_TYPE:-gene_id}"    # GTF attribute for grouping (column 9)
MIN_MAPQ="${MIN_MAPQ:-10}"                     # Minimum mapping quality filter

# ------------------------------------------------------------------------------
# 6. Logging Helpers
# ------------------------------------------------------------------------------
log_info()  { echo -e "\033[1;34m[INFO]\033[0m  $*"; }
log_step()  { echo -e "\033[1;35m[STEP]\033[0m  $*"; }
log_ok()    { echo -e "\033[1;32m[SUCCESS]\033[0m $*"; }
log_warn()  { echo -e "\033[1;33m[WARN]\033[0m  $*"; }
log_error() { echo -e "\033[1;31m[ERROR]\033[0m $*" >&2; }
