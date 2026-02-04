#!/usr/bin/env bash
set -euo pipefail

FILTERED_DIR=$1
NANOPLOT_DIR=$2

OUTDIR="MultiQC_report"

echo "=== MODULE 2: Running MultiQC ==="
echo "MultiQC output directory: ${OUTDIR}"

mkdir -p "${OUTDIR}"

multiqc \
    "${FILTERED_DIR}" \
    "${NANOPLOT_DIR}" \
    -o "${OUTDIR}" \
    --force
