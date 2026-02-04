#!/usr/bin/env python3
import os
import pandas as pd
from Bio import SeqIO
import subprocess
import gzip
import sys
import argparse

# ===============================
# 🔬 MÓDULO 3 – CONSENSO POR MUESTRA
# ===============================

parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True)
parser.add_argument('--fastq', required=True)
parser.add_argument('--assignments', required=True)
parser.add_argument('--racon_rounds', type=int, default=4)
args = parser.parse_args()

# Herramientas (Nextflow las buscará en el PATH del entorno conda)
VSEARCH = "vsearch"
MINIMAP2 = "minimap2"
RACON = "racon"
MIN_READS_FOR_RACON = 10

# 1. Cargar lecturas de la muestra en un diccionario
all_reads = {}
handle_open = gzip.open if args.fastq.endswith(".gz") else open
with handle_open(args.fastq, "rt") as f:
    for record in SeqIO.parse(f, "fastq"):
        all_reads[record.id] = record

# 2. Cargar asignaciones
df = pd.read_csv(args.assignments, sep="\t")
clusters = [c for c in df["cluster"].unique() if c != -1]

consensus_master = f"{args.sample}_clusters_consensus.fasta"

with open(consensus_master, "w") as master_handle:
    for cluster_id in clusters:
        cluster_reads_ids = df[df["cluster"] == cluster_id]["read_id"].tolist()
        if len(cluster_reads_ids) == 0: continue

        # Extraer reads del cluster
        tmp_fastq = f"tmp_c{cluster_id}.fastq"
        with open(tmp_fastq, "w") as f:
            for rid in cluster_reads_ids:
                if rid in all_reads:
                    SeqIO.write(all_reads[rid], f, "fastq")

        # 3.1 VSEARCH (Consenso inicial)
        tmp_vsearch = f"v_c{cluster_id}.fasta"
        subprocess.run([VSEARCH, "--cluster_fast", tmp_fastq, "--id", "0.51", "--consout", tmp_vsearch], 
                       capture_output=True, check=True)

        current_cons = tmp_vsearch
        
# 3.2 Racon Iterativo
        if len(cluster_reads_ids) >= MIN_READS_FOR_RACON:
            for i in range(1, args.racon_rounds + 1):
                tmp_sam = f"c{cluster_id}_r{i}.sam"
                tmp_racon = f"c{cluster_id}_r{i}.fasta"
                
                # VERIFICACIÓN DE SEGURIDAD: Si el consenso anterior está vacío, abortar pulido
                if not os.path.exists(current_cons) or os.path.getsize(current_cons) == 0:
                    print(f"⚠️ Alerta: Consenso previo vacío en cluster {cluster_id} ronda {i}. Usando último válido.")
                    break

                try:
                    with open(tmp_sam, "w") as sam:
                        subprocess.run([MINIMAP2, "-ax", "map-ont", current_cons, tmp_fastq], 
                                       stdout=sam, check=True, stderr=subprocess.DEVNULL)
                    
                    with open(tmp_racon, "w") as rac:
                        subprocess.run([RACON, tmp_fastq, tmp_sam, current_cons], 
                                       stdout=rac, check=True, stderr=subprocess.DEVNULL)
                    
                    # Si Racon generó un archivo válido, lo actualizamos
                    if os.path.exists(tmp_racon) and os.path.getsize(tmp_racon) > 0:
                        current_cons = tmp_racon
                    else:
                        break # Si Racon falla en generar algo, nos quedamos con la ronda anterior
                except subprocess.CalledProcessError:
                    print(f"⚠️ Error en pulido de cluster {cluster_id} ronda {i}. Saltando pulido extra.")
                    break
 
       # Guardar al maestro
        if os.path.exists(current_cons):
            for record in SeqIO.parse(current_cons, "fasta"):
                record.id = f"{args.sample}_c{cluster_id}"
                record.description = f"reads:{len(cluster_reads_ids)}"
                SeqIO.write(record, master_handle, "fasta")
        
        # Limpieza de archivos temporales del cluster
        for f in [tmp_fastq, tmp_vsearch] + [f"c{cluster_id}_r{i}.sam" for i in range(1, 5)]:
            if os.path.exists(f): os.remove(f)

print(f"✅ Consenso completado para {args.sample}")
