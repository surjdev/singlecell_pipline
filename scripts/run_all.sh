#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/run_all.sh
# Purpose: Master Runner for Production Pure-Bash Single-Cell Upstream Pipeline
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -i, --indir DIR        Directory with raw paired-end FASTQ reads (default: ${RAW_DIR})
  -o, --outdir DIR       Root output directory (default: ${PROJECT_ROOT})
  -g, --gtf PATH         Gene annotation GTF file (default: ${GTF_FILE})
  -x, --star-index DIR   STAR genome index directory (default: ${STAR_INDEX})
  -t, --threads NUM      CPU threads to use (default: ${THREADS})
  -c, --config FILE      Path to custom config.env file
  -h, --help             Show this help message and exit

Examples:
  # Run with defaults or config.env:
  bash $(basename "$0")

  # Run with custom parameters:
  bash $(basename "$0") -i /data/my_cells -x /refs/star_index -g /refs/genes.gtf -t 16
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--indir)      RAW_DIR="$2"; shift 2 ;;
        -o|--outdir)     PROJECT_ROOT="$2"; shift 2 ;;
        -g|--gtf)        GTF_FILE="$2"; shift 2 ;;
        -x|--star-index) STAR_INDEX="$2"; shift 2 ;;
        -t|--threads)    THREADS="$2"; shift 2 ;;
        -c|--config)     source "$2"; shift 2 ;;
        -h|--help)       usage ;;
        *) log_error "Unknown option: $1"; usage ;;
    esac
done

START_TIME=$(date +%s)

echo "=============================================================================="
echo " 🧬 PRODUCTION SMART-SEQ2 SINGLE-CELL UPSTREAM BASH PIPELINE"
echo "=============================================================================="
echo " Raw FASTQ directory : ${RAW_DIR}"
echo " Clean directory     : ${CLEAN_DIR}"
echo " Aligned directory   : ${ALIGNED_DIR}"
echo " Counts directory    : ${COUNTS_DIR}"
echo " Reports directory   : ${REPORTS_DIR}"
echo " STAR Index          : ${STAR_INDEX}"
echo " GTF Annotation      : ${GTF_FILE}"
echo " Execution Threads   : ${THREADS}"
echo " Start Timestamp     : $(date)"
echo "=============================================================================="

# 1. Verify required bioinformatics tools exist in PATH
REQUIRED_TOOLS=("fastqc" "fastp" "STAR" "samtools" "featureCounts" "multiqc")
MISSING_TOOLS=()
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        MISSING_TOOLS+=("${tool}")
    fi
done

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    log_error "Missing required bioinformatics tools: ${MISSING_TOOLS[*]}"
    log_info "Please ensure the pixi environment is activated or tools are installed:"
    log_info "  pixi shell"
    exit 1
fi

# 2. Stage 00: Raw FastQC
bash "${SCRIPT_DIR}/qc_data_00.sh" "${RAW_DIR}" "${QC_RAW_DIR}" "${THREADS}"

# 3. Stage 01: fastp Adapter & Quality Trimming
bash "${SCRIPT_DIR}/fastp_01.sh" "${RAW_DIR}" "${CLEAN_DIR}" "${FASTP_DIR}" "${THREADS}"

# 4. Stage 02: STAR Splice-Aware Alignment
bash "${SCRIPT_DIR}/align_02.sh" "${CLEAN_DIR}" "${ALIGNED_DIR}" "${STAR_INDEX}" "${THREADS}"

# 5. Stage 03: Subread featureCounts Quantification
bash "${SCRIPT_DIR}/quant_03.sh" "${ALIGNED_DIR}" "${COUNTS_DIR}" "${GTF_FILE}" "${THREADS}"

# 6. Stage 04: MultiQC Summary Aggregation
bash "${SCRIPT_DIR}/multiqc_04.sh" "${MULTIQC_DIR}" "${QC_RAW_DIR}" "${FASTP_DIR}" "${ALIGNED_DIR}" "${COUNTS_DIR}"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "=============================================================================="
echo " ✨ Pipeline Completed Successfully in ${ELAPSED} seconds!"
echo " Expression matrix: ${COUNTS_DIR}/gene_cell_count_matrix.tsv"
echo " MultiQC report:    ${MULTIQC_DIR}/single_cell_multiqc_report.html"
echo "=============================================================================="
