#!/usr/bin/env python3
import os
import gzip
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import Normalizer
from sklearn.decomposition import TruncatedSVD
import umap
import hdbscan
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# CONFIGURACIÓN Y ARGUMENTOS
# ===============================
parser = argparse.ArgumentParser(description='Module 2: Clustering with PCA+UMAP+HDBSCAN')
parser.add_argument('--input', required=True, help='Input FASTQ file (filtered)')
parser.add_argument('--sample', required=True, help='Sample name')
parser.add_argument('--min_len', type=int, default=200)
parser.add_argument('--max_len', type=int, default=400)
parser.add_argument('--pca_comp', type=int, default=500)
args = parser.parse_args()

# Variables globales que faltaban
KMER_SIZE = 7
MIN_CLUSTER_SIZE = 10

# ==============================
# 1️⃣ Leer y Filtrar secuencias
# ==============================
records = []
qc_removed = 0
print(f"=== Step 1: Loading reads for {args.sample} ===")

handle = gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "r")
for record in SeqIO.parse(handle, "fastq"):
    seq_len = len(record.seq)
    if args.min_len <= seq_len <= args.max_len:
        records.append((args.sample, record.id, str(record.seq).upper()))
    else:
        qc_removed += 1
handle.close()

if not records:
    print(f"⚠️ No reads found for sample {args.sample} after length filtering.")
    exit(0)

df = pd.DataFrame(records, columns=["sample", "read_id", "sequence"])
print(f"Loaded {len(df)} reads. Removed by QC: {qc_removed}")

# ==============================
# 2️⃣ Vectorización TF-IDF
# ==============================
print("=== Step 2: Computing TF-IDF weighted k-mer vectors ===")
def get_kmers(seq, k=KMER_SIZE):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

corpus = df["sequence"].tolist()
vectorizer = TfidfVectorizer(analyzer=lambda x: get_kmers(x, KMER_SIZE))
X = vectorizer.fit_transform(corpus)

normalizer = Normalizer(norm='l2')
X_norm = normalizer.fit_transform(X)

# ==============================
# 3️⃣ Reducción con PCA (SVD)
# ==============================
print(f"=== Step 3: Running PCA (TruncatedSVD) with {args.pca_comp} components ===")
n_comp = min(args.pca_comp, X_norm.shape[1] - 1)
svd = TruncatedSVD(n_components=n_comp, random_state=42)
X_pca = svd.fit_transform(X_norm)

# ==============================
# 4️⃣ UMAP
# ==============================
print("=== Step 4: Running UMAP ===")
umap_model = umap.UMAP(n_neighbors=20, min_dist=0.05, n_components=2, metric="cosine", random_state=42)
embedding = umap_model.fit_transform(X_pca)
df["UMAP1"] = embedding[:, 0]
df["UMAP2"] = embedding[:, 1]

# ==============================
# 5️⃣ HDBSCAN
# ==============================
print("=== Step 5: Clustering with HDBSCAN ===")
clusterer = hdbscan.HDBSCAN(min_cluster_size=MIN_CLUSTER_SIZE, min_samples=3, metric="euclidean")
df["cluster"] = clusterer.fit_predict(embedding)
n_clusters = len(set(df["cluster"])) - (1 if -1 in df["cluster"] else 0)

# ==============================
# 6️⃣ Guardar Resultados
# ==============================
df.to_csv(f"{args.sample}_read_cluster_assignments.tsv", sep="\t", index=False)

otu_table = df[df["cluster"] != -1].groupby(["cluster", "sample"]).size().unstack(fill_value=0)
otu_table.to_csv(f"{args.sample}_OTU_table.tsv", sep="\t")

summary = df["cluster"].value_counts().rename_axis("cluster").reset_index(name="num_reads")
summary.to_csv(f"{args.sample}_cluster_summary.tsv", sep="\t", index=False)

# ==============================
# 7️⃣ Gráfico de Eficiencia
# ==============================
noise_count = len(df[df["cluster"] == -1])
clustered_count = len(df[df["cluster"] != -1])
total = len(df)

data = pd.DataFrame({
    'Categoría': ['Agrupadas', 'Ruido'],
    'Porcentaje': [clustered_count/total*100, noise_count/total*100]
})

plt.figure(figsize=(8, 6))
sns.barplot(x='Categoría', y='Porcentaje', data=data)
plt.title(f'Eficiencia Clustering: {args.sample}')
plt.savefig(f"{args.sample}_clustering_efficiency_plot.png")
plt.close()

print(f"=== ✅ Finished {args.sample}. Clusters found: {n_clusters} ===")
