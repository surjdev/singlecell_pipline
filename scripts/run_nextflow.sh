#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/run_nextflow.sh
# Purpose: Direct Runner for Single-Cell Nextflow Pipeline (H100 / GPU / Dedicated Server)
# ==============================================================================
# 🏷️ Usage:
#   bash scripts/run_nextflow.sh [OPTIONS]
#
# Options:
#   --threads, -t <num> : Number of CPU threads (default: auto-detected from nproc)
#   --conda, -c <path>  : Path to custom conda environment (optional)
#   --bg                : Run in background using nohup (safe from terminal disconnect)
#   --help, -h          : Show this help message
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

# Default parameters
THREADS="$(nproc 2>/dev/null || echo 16)"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-}"
RUN_BG=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--conda)
            CONDA_ENV_PATH="$2"
            shift 2
            ;;
        --bg)
            RUN_BG=true
            shift
            ;;
        -h|--help)
            echo "Usage: bash scripts/run_nextflow.sh [--threads <N>] [--conda <path>] [--bg]"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1"
            exit 1
            ;;
    esac
done

# If requested to run in background, re-launch via nohup
if [ "$RUN_BG" = true ]; then
    mkdir -p logs
    LOG_FILE="logs/nextflow_$(date +%Y%m%d_%H%M%S).log"
    echo "======================================================================"
    echo " Launching Nextflow in BACKGROUND (Safe from browser/terminal closure)"
    echo " Log file: ${LOG_FILE}"
    echo " To view progress: tail -f ${LOG_FILE}"
    echo " To check process: ps aux | grep nextflow"
    echo "======================================================================"
    nohup bash "$0" --threads "${THREADS}" ${CONDA_ENV_PATH:+--conda "${CONDA_ENV_PATH}"} > "${LOG_FILE}" 2>&1 &
    exit 0
fi

echo "======================================================================"
echo " Starting Smart-seq2 Nextflow Pipeline"
echo " Host Machine:    $(hostname)"
echo " CPU Threads:     ${THREADS}"
echo " Working Dir:     ${PROJECT_ROOT}"
echo " Start Time:      $(date)"
echo "======================================================================"

mkdir -p logs work reports

# 1. Activate Conda Environment (if specified or found)
if [ -n "${CONDA_ENV_PATH}" ] && [ -d "${CONDA_ENV_PATH}" ]; then
    if [ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]; then
        # shellcheck source=/dev/null
        source "$(conda info --base)/etc/profile.d/conda.sh"
    fi
    echo "[INFO] Activating Conda: ${CONDA_ENV_PATH}"
    conda activate "${CONDA_ENV_PATH}"
fi

# 2. Check Nextflow availability
if ! command -v nextflow >/dev/null 2>&1; then
    echo "[ERROR] 'nextflow' is not in PATH! Please activate your conda environment first."
    exit 1
fi

# 3. Synchronize bioinformatics dependencies via Pixi (if present)
if command -v pixi >/dev/null 2>&1 && [ -f "pixi.toml" ]; then
    echo "[INFO] Syncing tools with pixi..."
    pixi install
fi

# 4. Storage & Cache Redirection
export NXF_WORK="${PROJECT_ROOT}/work"
export NXF_TEMP="${PROJECT_ROOT}/work/tmp"
export TMPDIR="${PROJECT_ROOT}/work/tmp"
mkdir -p "${NXF_WORK}" "${NXF_TEMP}"

# 5. Execute Nextflow Pipeline
echo "[STEP] Executing Nextflow..."
nextflow run nextflow/main.nf \
    -profile local \
    --threads "${THREADS}" \
    -resume

echo "======================================================================"
echo " Pipeline Finished Successfully!"
echo " End Time:       $(date)"
echo " Results:        ${PROJECT_ROOT}/results"
echo " MultiQC Report: ${PROJECT_ROOT}/results/multiqc/multiqc_report.html"
echo "======================================================================"
