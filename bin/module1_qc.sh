#!/usr/bin/env bash
set -euo pipefail

# =========================
# Parse arguments
# =========================
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample) SAMPLE="$2"; shift 2 ;;
        --input_dir) INPUT_DIR="$2"; shift 2 ;; # Nextflow pasará el path de la carpeta
        --primer_f) PRIMER_F="$2"; shift 2 ;;
        --primer_r) PRIMER_R="$2"; shift 2 ;;
        *) echo "❌ Unknown option $1"; exit 1 ;;
    esac
done

echo "============================================"
echo "🧬 MODULE 1: QC & CONCATENATION"
echo "🧪 Sample   : $SAMPLE"
echo "📂 Input    : $INPUT_DIR"
echo "============================================"

# ==================================================
# 1. Concatenación con FASTCAT (Requisito Especial)
# ==================================================
# Fastcat es inteligente: maneja .gz y .fastq mezclados 
# y genera un reporte de estadísticas inicial.
OUT_CONCAT="${SAMPLE}_concatenated.fastq.gz"

echo "📥 Running fastcat..."
fastcat "$INPUT_DIR" | gzip -c > "$OUT_CONCAT"

# ==================================================
# 2. Primer trimming con CUTADAPT
# ==================================================
echo "✂️ Running cutadapt..."
cutadapt \
    -g "$PRIMER_F" \
    -a "$PRIMER_R" \
    --revcomp \
    -e 0.20 \
    --discard-untrimmed \
    --minimum-length 200 \
    --maximum-length 500 \
    -o "${SAMPLE}_V4_filtered.fastq.gz" \
    "$OUT_CONCAT" \
    > "${SAMPLE}_cutadapt.log" 2>&1 || true

# Garantizar que el archivo existe aunque esté vacío (evita que NF falle)
if [[ ! -s "${SAMPLE}_V4_filtered.fastq.gz" ]]; then
    echo "⚠️ No reads passed filtering"
    echo "" | gzip -c > "${SAMPLE}_V4_filtered.fastq.gz"
fi

# ==================================================
# 3. Reporte de Calidad con NANOPLOT
# ==================================================
if command -v NanoPlot >/dev/null 2>&1; then
    echo "📊 Running NanoPlot..."
    NanoPlot \
        --fastq "${SAMPLE}_V4_filtered.fastq.gz" \
        --outdir "nanoplot_${SAMPLE}" \
        --prefix "${SAMPLE}_" \
        --threads 4 \
        --plots hex
fi

echo "✅ QC DONE: $SAMPLE"
