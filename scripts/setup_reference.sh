#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/setup_reference.sh
# Purpose: Download Standard Reference Genomes and Build Production STAR Index
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  --species SPECIES      Species to download: 'human' (GRCh38) or 'mouse' (GRCm39)
  --fasta PATH           Path to custom genome FASTA (skip download)
  --gtf PATH             Path to custom gene annotation GTF (skip download)
  --outdir DIR           Output directory for reference and index (default: data/reference)
  --threads NUM          Number of CPU threads for STAR indexing (default: ${THREADS})
  --sjdbOverhang NUM     Read length - 1 for splice junctions (default: 100)
  -h, --help             Show this help message and exit

Examples:
  # Download human GENCODE GRCh38 and build index:
  bash $(basename "$0") --species human --threads 16

  # Build STAR index from existing local FASTA and GTF:
  bash $(basename "$0") --fasta /path/to/genome.fa --gtf /path/to/genes.gtf --threads 16
EOF
    exit 0
}

SPECIES=""
CUSTOM_FASTA=""
CUSTOM_GTF=""
TARGET_DIR="${REF_DIR}"
SJDB_OVERHANG=100

while [[ $# -gt 0 ]]; do
    case "$1" in
        --species) SPECIES="$2"; shift 2 ;;
        --fasta)   CUSTOM_FASTA="$2"; shift 2 ;;
        --gtf)     CUSTOM_GTF="$2"; shift 2 ;;
        --outdir)  TARGET_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --sjdbOverhang) SJDB_OVERHANG="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) log_error "Unknown option: $1"; usage ;;
    esac
done

mkdir -p "${TARGET_DIR}"
TARGET_INDEX="${TARGET_DIR}/star_index"

# 1. Obtain Genome FASTA and GTF
if [ -n "${CUSTOM_FASTA}" ] && [ -n "${CUSTOM_GTF}" ]; then
    log_info "Using user-provided FASTA: ${CUSTOM_FASTA}"
    log_info "Using user-provided GTF:   ${CUSTOM_GTF}"
    FASTA_FILE="${CUSTOM_FASTA}"
    GTF_INPUT="${CUSTOM_GTF}"
elif [ "${SPECIES}" == "human" ]; then
    log_step "Downloading Human (GRCh38 / GENCODE release 44) reference..."
    FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
    GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
    
    FASTA_GZ="${TARGET_DIR}/GRCh38.primary_assembly.genome.fa.gz"
    GTF_GZ="${TARGET_DIR}/gencode.v44.annotation.gtf.gz"
    FASTA_FILE="${TARGET_DIR}/genome.fa"
    GTF_INPUT="${TARGET_DIR}/genes.gtf"

    if [ ! -f "${FASTA_FILE}" ]; then
        log_info "Downloading genome FASTA from ${FASTA_URL}..."
        curl -L -o "${FASTA_GZ}" "${FASTA_URL}"
        log_info "Decompressing genome FASTA..."
        gzip -dc "${FASTA_GZ}" > "${FASTA_FILE}"
        rm -f "${FASTA_GZ}"
    fi

    if [ ! -f "${GTF_INPUT}" ]; then
        log_info "Downloading gene GTF from ${GTF_URL}..."
        curl -L -o "${GTF_GZ}" "${GTF_URL}"
        log_info "Decompressing gene GTF..."
        gzip -dc "${GTF_GZ}" > "${GTF_INPUT}"
        rm -f "${GTF_GZ}"
    fi
elif [ "${SPECIES}" == "mouse" ]; then
    log_step "Downloading Mouse (GRCm39 / GENCODE release M33) reference..."
    FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz"
    GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz"

    FASTA_GZ="${TARGET_DIR}/GRCm39.primary_assembly.genome.fa.gz"
    GTF_GZ="${TARGET_DIR}/gencode.vM33.annotation.gtf.gz"
    FASTA_FILE="${TARGET_DIR}/genome.fa"
    GTF_INPUT="${TARGET_DIR}/genes.gtf"

    if [ ! -f "${FASTA_FILE}" ]; then
        log_info "Downloading genome FASTA from ${FASTA_URL}..."
        curl -L -o "${FASTA_GZ}" "${FASTA_URL}"
        log_info "Decompressing genome FASTA..."
        gzip -dc "${FASTA_GZ}" > "${FASTA_FILE}"
        rm -f "${FASTA_GZ}"
    fi

    if [ ! -f "${GTF_INPUT}" ]; then
        log_info "Downloading gene GTF from ${GTF_URL}..."
        curl -L -o "${GTF_GZ}" "${GTF_URL}"
        log_info "Decompressing gene GTF..."
        gzip -dc "${GTF_GZ}" > "${GTF_INPUT}"
        rm -f "${GTF_GZ}"
    fi
else
    log_error "Please specify either --species (human|mouse) or both --fasta and --gtf!"
    usage
fi

# 2. Build Production STAR Index
mkdir -p "${TARGET_INDEX}"
log_step "Building Production STAR Index in ${TARGET_INDEX}..."
log_info "Threads: ${THREADS} | sjdbOverhang: ${SJDB_OVERHANG}"

STAR \
    --runMode genomeGenerate \
    --genomeDir "${TARGET_INDEX}" \
    --genomeFastaFiles "${FASTA_FILE}" \
    --sjdbGTFfile "${GTF_INPUT}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    --runThreadN "${THREADS}"

log_ok "Reference genome and STAR index successfully built at ${TARGET_INDEX}!"
EOF
