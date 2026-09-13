#!/usr/bin/env bash
# ==============================================================================
# SLURM Submission Script for Smart-seq2 Single-Cell Nextflow Pipeline
# ==============================================================================
#SBATCH --job-name=sc_smartseq2
#SBATCH --partition=compute              # Adjust to your cluster partition (e.g. batch, normal, compute)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16               # CPUs allocated for multithreading (STAR, fastp, featureCounts)
#SBATCH --mem=64G                        # >= 50-64GB required for Human GRCh38 STAR alignment
#SBATCH --time=12:00:00                  # Maximum run time
#SBATCH --output=logs/slurm_%j.out       # Standard output log
#SBATCH --error=logs/slurm_%j.err        # Standard error log

set -euo pipefail

echo "======================================================================"
echo " Starting Single-Cell Pipeline on HPC (SLURM Job ID: ${SLURM_JOB_ID:-Manual})"
echo " Running on Node: $(hostname)"
echo " Start Time:      $(date)"
echo " Cores Allocated: ${SLURM_CPUS_PER_TASK:-16}"
echo " Memory:          64GB"
echo "======================================================================"

# 1. Project directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"
mkdir -p logs "${PROJECT_ROOT}/work" "${PROJECT_ROOT}/reports"

# 2. Activate Conda Environment (Adjust CONDA_ENV_PATH to your install prefix)
CONDA_ENV_PATH="${CONDA_ENV_PATH:-/home/koraop/nob_dir/conda_envs/sc_env}"

if [ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi

if [ -d "${CONDA_ENV_PATH}" ]; then
    echo "[INFO] Activating Conda environment from: ${CONDA_ENV_PATH}"
    conda activate "${CONDA_ENV_PATH}"
else
    echo "[WARN] Conda prefix ${CONDA_ENV_PATH} not found. Trying active environment..."
fi

# Verify core orchestration tools
echo "[INFO] Verifying tools:"
command -v nextflow >/dev/null 2>&1 && nextflow -version || { echo "[ERROR] nextflow not found!"; exit 1; }
command -v pixi >/dev/null 2>&1 && pixi --version || echo "[WARN] pixi not in PATH, checking conda env tools..."

# 3. Synchronize bioinformatics dependencies via Pixi (if using pixi environment)
if command -v pixi >/dev/null 2>&1 && [ -f "pixi.toml" ]; then
    echo "[INFO] Syncing bioinformatics tools via pixi..."
    pixi install
fi

# 4. Storage & Cache Redirection for HPC (Avoid filling up /tmp and $HOME)
export NXF_WORK="${PROJECT_ROOT}/work"
export NXF_TEMP="${PROJECT_ROOT}/work/tmp"
export TMPDIR="${PROJECT_ROOT}/work/tmp"
mkdir -p "${NXF_WORK}" "${NXF_TEMP}"

# 5. Run Nextflow Pipeline with Resume
# Profile 'local' runs within this allocated 16-core / 64GB compute node efficiently
echo "[STEP] Launching Nextflow Pipeline..."
nextflow run nextflow/main.nf \
    -profile local \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    -resume

echo "======================================================================"
echo " Single-Cell Pipeline Completed Successfully!"
echo " Finished Time: $(date)"
echo " Results available in: ${PROJECT_ROOT}/results"
echo " MultiQC report:       ${PROJECT_ROOT}/results/multiqc/multiqc_report.html"
echo "======================================================================"
