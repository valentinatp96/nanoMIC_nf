#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# ===============================
# CONFIGURACIÓN DE ARGUMENTOS (Nextflow style)
# ===============================
if len(sys.argv) != 3:
    print("Uso: python3 module5_visualizations.py <OTU_TABLE_TSV> <OUTPUT_DIR>")
    sys.exit(1)

FINAL_OTU_TSV = sys.argv[1]
PLOTS_DIR = sys.argv[2]
os.makedirs(PLOTS_DIR, exist_ok=True)

# ==============================
# CARGAR DATOS
# ==============================
try:
    otu_df = pd.read_csv(FINAL_OTU_TSV, sep="\t")
    # Si 'cluster' es una columna, la ponemos como índice
    if 'cluster' in otu_df.columns:
        otu_df.set_index('cluster', inplace=True)
except Exception as e:
    print(f"❌ ERROR: No se pudo cargar la tabla: {e}")
    sys.exit(1)

# Identificar muestras (todas menos la columna Taxonomy)
samples = [c for c in otu_df.columns if c != "Taxonomy"]
tax_col = "Taxonomy"

# ==============================
# 1️⃣ Gráficas de taxonomía
# ==============================
tax_series = otu_df[tax_col].fillna('Unassigned').copy()
tax_df = tax_series.str.split('; ', expand=True)

tax_levels_full = ['Domain','Phylum','Class','Order','Family','Genus','Species']
tax_df.columns = tax_levels_full[:tax_df.shape[1]]

def plot_tax_level(level, otu_data, tax_data, outdir):
    tax_level_series = tax_data[level].str.replace(f'{level}__', '', regex=False).fillna(f'Unassigned_{level}')
    
    df_processing = otu_data[samples].copy()
    df_processing['Taxon_Level'] = tax_level_series.values
    df_sum = df_processing.groupby('Taxon_Level')[samples].sum()
    
    # Normalizar a abundancia relativa
    df_rel = df_sum.div(df_sum.sum(axis=0), axis=1) * 100

    top_n = 20
    top_taxa = df_rel.sum(axis=1).nlargest(top_n).index
    df_top = df_rel.loc[top_taxa]
    df_other = df_rel.drop(top_taxa, errors='ignore').sum(axis=0)
    df_final = pd.concat([df_top, pd.DataFrame([df_other], index=['Others'])])

    plt.figure(figsize=(14,10))
    df_final.T.plot(kind='bar', stacked=True, width=0.8, colormap='tab20')
    plt.ylabel('Abundancia relativa (%)')
    plt.title(f'Composición taxonómica Nivel: {level}')
    plt.legend(bbox_to_anchor=(1.05,1), loc='upper left', fontsize='small')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{level}_abundance.png"), dpi=300)
    plt.close()

levels_to_plot = [lvl for lvl in ['Phylum','Class','Order','Family','Genus'] if lvl in tax_df.columns]
for lvl in levels_to_plot:
    plot_tax_level(lvl, otu_df, tax_df, PLOTS_DIR)

# ==============================
# 2️⃣ Curvas de rarefacción
# ==============================
def rarefaction_curve(counts, step=100, iterations=10):
    counts = np.array(counts)
    population = np.repeat(np.arange(len(counts)), counts)
    total_reads = len(population)
    
    if total_reads < step: return [0], [0]
    
    points = np.arange(step, total_reads, step)
    results = []
    for _ in range(iterations):
        obs = [len(np.unique(np.random.choice(population, n, replace=False))) for n in points]
        results.append(obs)
    return points, np.mean(results, axis=0)

plt.figure(figsize=(10,6))
for sample in samples:
    x, y = rarefaction_curve(otu_df[sample].values)
    plt.plot(x, y, label=sample)
plt.xlabel("Secuencias")
plt.ylabel("OTUs observados")
plt.legend(bbox_to_anchor=(1.05,1))
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "rarefaction_curves.png"), dpi=300)
print(f"✅ Visualizaciones completadas en {PLOTS_DIR}")
