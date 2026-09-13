#!/usr/bin/env bash
# ==============================================================================
# Script: scripts/download_data.sh
# Purpose: Download Single-Cell FASTQ files from CSV with Resume capability
# ==============================================================================
# 🏷️ Features & Options:
#   --csv <file>     : Path to metadata CSV (default: single_cell_cls_balance.csv)
#   --outdir <dir>   : Directory to save downloaded FASTQ files (default: data/raw)
#   --link-sample    : Create symlink with biological sample name (e.g. HP1502401_A14.fastq.gz)
#   --ftp            : Force FTP protocol instead of converting to HTTPS (default: HTTPS)
#   --dry-run        : Preview files to be downloaded without actually fetching them
#   --help, -h       : Show this help message
#
# 🔁 Resume Capability:
#   - Uses 'wget -c' or 'curl -C -' to resume interrupted transfers.
#   - Automatically tests completed gzip files with 'gzip -t' to verify integrity.
#   - If file exists and passes integrity check, download is safely skipped.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Load environment configuration if available
if [ -f "${SCRIPT_DIR}/env.sh" ]; then
    # shellcheck source=/dev/null
    source "${SCRIPT_DIR}/env.sh"
else
    log_info()  { echo -e "\033[1;34m[INFO]\033[0m  $*"; }
    log_step()  { echo -e "\033[1;35m[STEP]\033[0m  $*"; }
    log_ok()    { echo -e "\033[1;32m[SUCCESS]\033[0m $*"; }
    log_warn()  { echo -e "\033[1;33m[WARN]\033[0m  $*"; }
    log_error() { echo -e "\033[1;31m[ERROR]\033[0m $*" >&2; }
    RAW_DIR="${PROJECT_ROOT}/data/raw"
fi

# Defaults
CSV_FILE="${PROJECT_ROOT}/single_cell_cls_balance.csv"
OUT_DIR="${RAW_DIR:-${PROJECT_ROOT}/data/raw}"
USE_HTTPS=true
LINK_SAMPLE=true
DRY_RUN=false

show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Download Single-Cell FASTQ files with automatic resume support.

Options:
  --csv <path>       Metadata CSV file (default: ${CSV_FILE})
  --outdir <path>    Output directory for FASTQ files (default: ${OUT_DIR})
  --no-symlink       Do not create sample-name symlinks (only keep ERR run IDs)
  --ftp              Keep ftp:// URLs (do not convert to https://)
  --dry-run          Print URLs without downloading
  -h, --help         Show this help message

Examples:
  # Standard execution with resume:
  bash scripts/download_data.sh

  # Custom CSV and output folder:
  bash scripts/download_data.sh --csv single_cell_cls_balance.csv --outdir data/raw
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --csv)
            CSV_FILE="$2"
            shift 2
            ;;
        --outdir)
            OUT_DIR="$2"
            shift 2
            ;;
        --no-symlink)
            LINK_SAMPLE=false
            shift
            ;;
        --ftp)
            USE_HTTPS=false
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validation
if [ ! -f "${CSV_FILE}" ]; then
    log_error "Metadata CSV not found: ${CSV_FILE}"
    exit 1
fi

mkdir -p "${OUT_DIR}"

log_step "Initializing FASTQ Downloader"
log_info "Metadata CSV:     ${CSV_FILE}"
log_info "Destination Dir:  ${OUT_DIR}"
log_info "Protocol:         $([ "$USE_HTTPS" = true ] && echo 'HTTPS (recommended/firewall safe)' || echo 'FTP')"
log_info "Create Symlink:   ${LINK_SAMPLE}"

# Detect download tool (wget or curl)
DOWNLOADER=""
if command -v wget >/dev/null 2>&1; then
    DOWNLOADER="wget"
elif command -v curl >/dev/null 2>&1; then
    DOWNLOADER="curl"
else
    log_error "Neither 'wget' nor 'curl' is installed. Please install one of them."
    exit 1
fi
log_info "Downloader Engine: ${DOWNLOADER} (with resume enabled)"

# Extract download entries (run_id, sample_name, submitted_name, fastq_uri)
# Preference: Use python3 for robust CSV parsing; fallback to awk
ENTRIES=()
if command -v python3 >/dev/null 2>&1; then
    while IFS=$'\t' read -r run_id sample_name sub_name uri; do
        if [ -n "$uri" ] && [ "$uri" != "Comment[FASTQ_URI]" ]; then
            ENTRIES+=("${run_id}	${sample_name}	${sub_name}	${uri}")
        fi
    done < <(python3 -c "
import csv, sys
with open(sys.argv[1], 'r', encoding='utf-8', errors='replace') as f:
    reader = csv.DictReader(f)
    for row in reader:
        run = row.get('Comment[ENA_RUN]', '').strip()
        sample = row.get('Source Name', '').strip()
        sub = row.get('Comment[SUBMITTED_FILE_NAME]', '').strip()
        uri = row.get('Comment[FASTQ_URI]', '').strip()
        if uri:
            print(f'{run}\t{sample}\t{sub}\t{uri}')
" "${CSV_FILE}")
else
    # Fallback to awk (col 44=ENA_RUN, col 2=Source Name, col 43=SUBMITTED_FILE_NAME, col 45=FASTQ_URI)
    while IFS=$'\t' read -r run_id sample_name sub_name uri; do
        if [ -n "$uri" ] && [ "$uri" != "Comment[FASTQ_URI]" ]; then
            ENTRIES+=("${run_id}	${sample_name}	${sub_name}	${uri}")
        fi
    done < <(awk -F',' 'NR > 1 { gsub(/\"/, "", $0); print $44 "\t" $2 "\t" $43 "\t" $45 }' "${CSV_FILE}")
fi

TOTAL_COUNT="${#ENTRIES[@]}"
if [ "$TOTAL_COUNT" -eq 0 ]; then
    log_warn "No valid FASTQ download URIs found in ${CSV_FILE}"
    exit 0
fi

log_info "Found ${TOTAL_COUNT} files to process."

# Counters
DOWNLOADED=0
SKIPPED=0
FAILED=0
CURRENT=0

for entry in "${ENTRIES[@]}"; do
    CURRENT=$((CURRENT + 1))
    IFS=$'\t' read -r run_id sample_name sub_name uri <<< "${entry}"

    # Target URL
    DOWNLOAD_URL="${uri}"
    if [ "$USE_HTTPS" = true ]; then
        DOWNLOAD_URL="${DOWNLOAD_URL/ftp:\/\//https:\/\/}"
    fi

    # Target file name from URL (e.g. ERR1630113.fastq.gz)
    FILE_NAME="$(basename "${DOWNLOAD_URL}")"
    TARGET_PATH="${OUT_DIR}/${FILE_NAME}"

    echo "----------------------------------------------------------------------"
    log_step "[${CURRENT}/${TOTAL_COUNT}] Sample: ${sample_name} (${run_id})"
    log_info "URL:  ${DOWNLOAD_URL}"
    log_info "File: ${TARGET_PATH}"

    if [ "$DRY_RUN" = true ]; then
        log_info "[DRY-RUN] Would download to: ${TARGET_PATH}"
        continue
    fi

    # Check if file already exists and is a valid non-empty gzip
    if [ -f "${TARGET_PATH}" ] && [ -s "${TARGET_PATH}" ]; then
        # Quick integrity test using gzip
        if gzip -t "${TARGET_PATH}" >/dev/null 2>&1; then
            log_ok "File already exists and passed integrity test (gzip -t). Skipping."
            SKIPPED=$((SKIPPED + 1))
            
            # Create symlink if requested
            if [ "$LINK_SAMPLE" = true ] && [ -n "${sub_name}" ]; then
                SAMPLE_LINK="${OUT_DIR}/${sub_name}"
                if [ ! -e "${SAMPLE_LINK}" ]; then
                    ln -sf "${FILE_NAME}" "${SAMPLE_LINK}"
                    log_info "Linked sample name: ${sub_name} -> ${FILE_NAME}"
                fi
            fi
            continue
        else
            log_warn "Corrupted or partial file detected. Will attempt to resume download..."
        fi
    fi

    # Perform download with resume support
    DOWNLOAD_SUCCESS=false
    if [ "$DOWNLOADER" = "wget" ]; then
        # wget -c: continue getting a partially-downloaded file
        if wget -c \
            --tries=10 \
            --timeout=30 \
            --waitretry=2 \
            --retry-connrefused \
            --show-progress \
            -O "${TARGET_PATH}" \
            "${DOWNLOAD_URL}"; then
            DOWNLOAD_SUCCESS=true
        fi
    elif [ "$DOWNLOADER" = "curl" ]; then
        # curl -C -: continue at offset from existing file
        if curl -C - \
            -L \
            --retry 5 \
            --retry-delay 2 \
            --retry-connrefused \
            --connect-timeout 30 \
            --fail \
            -o "${TARGET_PATH}" \
            "${DOWNLOAD_URL}"; then
            DOWNLOAD_SUCCESS=true
        fi
    fi

    # Verify download result
    if [ "$DOWNLOAD_SUCCESS" = true ] && [ -f "${TARGET_PATH}" ] && [ -s "${TARGET_PATH}" ]; then
        if gzip -t "${TARGET_PATH}" >/dev/null 2>&1; then
            log_ok "Download completed and verified successfully: ${FILE_NAME}"
            DOWNLOADED=$((DOWNLOADED + 1))

            # Symlink sample name (e.g. HP1502401_A14.fastq.gz -> ERR1630113.fastq.gz)
            if [ "$LINK_SAMPLE" = true ] && [ -n "${sub_name}" ]; then
                SAMPLE_LINK="${OUT_DIR}/${sub_name}"
                ln -sf "${FILE_NAME}" "${SAMPLE_LINK}"
                log_info "Linked sample name: ${sub_name} -> ${FILE_NAME}"
            fi
        else
            log_error "File was downloaded but failed integrity check (corrupted gzip): ${FILE_NAME}"
            FAILED=$((FAILED + 1))
        fi
    else
        log_error "Failed to download: ${DOWNLOAD_URL}"
        FAILED=$((FAILED + 1))
    fi
done

echo "======================================================================"
log_step "Download Summary"
log_info "Total Files:     ${TOTAL_COUNT}"
log_ok   "Downloaded:      ${DOWNLOADED}"
log_info "Skipped/Exists:  ${SKIPPED}"
if [ "$FAILED" -gt 0 ]; then
    log_error "Failed:          ${FAILED} (Re-run script to resume failed files)"
    exit 1
else
    log_ok "All files are ready in ${OUT_DIR}"
fi
