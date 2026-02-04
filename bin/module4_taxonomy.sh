#!/bin/bash
set -euo pipefail

# === 1. Captura de Argumentos (Recibidos desde el .nf) ===
RAW_FASTA=$1
ASSIGNMENTS=$2
DB_SEQS=$3
DB_TAX=$4

# Variables de trabajo locales (Nextflow trabaja en la carpeta actual)
RENAMED_FASTA="consensus_renamed.fasta"
CONSENSUS_QZA="consensus_renamed.qza"
TAXONOMY_QZA="taxonomy_vsearch.qza"
TAXONOMY_TSV="taxonomy_vsearch.tsv"
OUTPUT_FINAL="OTU_table_taxonomy.tsv"
QC_PLOT="taxonomy_assignment_qc.png"

echo "===================================================="
echo "🧬 Iniciando Asignación Taxonómica con QIIME2"
echo "===================================================="

# === STEP 1: Renombrar cabeceras del FASTA ===
# Esto evita errores en QIIME2 si hay nombres duplicados entre muestras
awk '/^>/{print ">cluster_" ++i; next}{print}' "${RAW_FASTA}" > "${RENAMED_FASTA}"

# === STEP 2: Importar a QIIME2 ===
qiime tools import \
    --input-path "${RENAMED_FASTA}" \
    --output-path "${CONSENSUS_QZA}" \
    --type 'FeatureData[Sequence]'

# === STEP 3: Clasificación VSEARCH ===
qiime feature-classifier classify-consensus-vsearch \
    --i-query "${CONSENSUS_QZA}" \
    --i-reference-reads "${DB_SEQS}" \
    --i-reference-taxonomy "${DB_TAX}" \
    --p-min-consensus 0.51 \
    --o-classification "${TAXONOMY_QZA}" \
    --o-search-results "search_results.qza"

# Exportar el TSV de taxonomía
qiime tools export --input-path "${TAXONOMY_QZA}" --output-path .
mv taxonomy.tsv "${TAXONOMY_TSV}"

# === STEP 4: Generar Tabla OTU Final (Python) ===
# Pasamos las rutas al bloque Python mediante variables de entorno
export ASSIGNMENTS
export TAXONOMY_TSV
export OUTPUT_FINAL
export QC_PLOT

python3 <<PYTHONCODE
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# A. Cargar datos
df_clusters = pd.read_csv(os.environ['ASSIGNMENTS'], sep="\t")
df_clusters = df_clusters[df_clusters['cluster'] != -1] # Quitar ruido

# B. Crear tabla de abundancia (OTU Table)
otu = df_clusters.groupby(["cluster", "sample"]).size().unstack(fill_value=0)

# C. Mapear nombres de clusters para que coincidan con el FASTA (cluster_1, 2...)
# Importante: El orden del groupby debe coincidir con el orden del awk anterior
cluster_id_map = {id: f"cluster_{i+1}" for i, id in enumerate(otu.index)}
otu.index = otu.index.map(cluster_id_map)
otu.index.name = 'cluster'

# D. Cargar Taxonomía
tax = pd.read_csv(os.environ['TAXONOMY_TSV'], sep="\t")
tax = tax.rename(columns={'Feature ID':'cluster', 'Taxon':'Taxonomy'})
tax['cluster'] = tax['cluster'].astype(str)
otu.index = otu.index.astype(str)

# E. Merge y Guardar
merged = otu.merge(tax[['cluster','Taxonomy']], on='cluster', how='left')
merged.to_csv(os.environ['OUTPUT_FINAL'], sep="\t", index=False)

# F. Generar QC Plot
assigned = tax['Taxonomy'].str.contains(';', na=False).sum()
unassigned = len(tax) - assigned
plt.figure(figsize=(6,4))
sns.barplot(x=['Asignado', 'No Asignado'], y=[assigned, unassigned], palette='viridis')
plt.title(f'Calidad de Asignación (Total: {len(tax)} OTUs)')
plt.savefig(os.environ['QC_PLOT'])
PYTHONCODE

echo "✅ Módulo 4 finalizado con éxito."
