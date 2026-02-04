# nanoMIC-nf

### Workflow de Nextflow para la caracterización de microbiomas asociados a corrosión

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Resumen Técnico
**nanoMIC-nf** es un pipeline bioinformático desarrollado en **Nextflow** diseñado para procesar datos de secuenciación de tercera generación (Oxford Nanopore Technologies). El flujo de trabajo implementa una arquitectura modular para la caracterización taxonómica de comunidades microbianas implicadas en procesos de biocorrosión, garantizando la reproducibilidad mediante el uso de contenedores y una gestión eficiente de recursos.

Este proyecto constituye el núcleo técnico de mi trabajo de fin de máster (TFM).

---

## Estructura del Repositorio
El repositorio sigue las mejores prácticas de desarrollo de Nextflow (DSL2):

* `main.nf`: Punto de entrada principal para la ejecución del flujo.
* `workflows/`: Lógica de orquestación del pipeline (`nanopore_workflow.nf`).
* `modules/`: Definición de procesos independientes para cada etapa del análisis.
* `bin/`: Scripts personalizados (Bash y Python) que ejecutan la lógica técnica.
* `nextflow.config`: Archivo de configuración global y gestión de perfiles.

---

## Arquitectura de Módulos

### 1. Control de Calidad (QC)
* **Herramientas:** `Fastcat`, `Cutadapt` `NanoPlot`.
* **Función:** Concatenación de archivos FASTQ crudos, generación de estadísticas de calidad (longitud de lectura, calidad base N50) y filtrado por umbrales de longitud.

### 2. Clustering 
* **Herramientas:** Scripts personalizados en Python (`module2_clustering.py`).
* **Función:** Agrupamiento de lecturas en Unidades Taxonómicas Operativas (OTUs) y discriminación de ruido técnico para limpiar la señal biológica.

### 3. Generación de Consenso
* **Función:** Cálculo de secuencias consenso por cada cluster para mitigar la tasa de error característica de Nanopore, mejorando la fidelidad de la secuencia final.

### 4. Asignación Taxonómica
* **Herramientas:** `QIIME2`, `VSEARCH`.
* **Función:** Integración de secuencias consenso en el ecosistema QIIME2. Utiliza la base de datos **SILVA 138** para la clasificación taxonómica mediante alineamiento y consenso de identidad.

### 5. Visualización y Análisis de Biodiversidad
* **Herramientas:** Python (`Pandas`, `Seaborn`).
* **Función:** Automatización de curvas de rarefacción y gráficos de abundancia relativa desde el nivel de Phylum hasta Género.

### 6. Reporte Global (MultiQC)
* **Función:** Agregación de métricas de todos los módulos en un reporte interactivo final para validación rápida.

---

## Guía de Uso

### Requisitos previos
- Nextflow (v22.10.0 o superior)
- Docker, Singularity o Conda (recomendado para reproducibilidad)

### Ejecución básica
```bash
nextflow run main.nf \
  --input 'ruta/a/Input-Data/raw_data/' \
  --databases 'ruta/a/qiime2_databases/' \
  --outdir 'results' \
  -profile docker

---

## Información Académica
**Tesis:** Automatización de un flujo de trabajo para la caracterización de microbiomas asociados a corrosión  
**Autora:** Valentina Tapia Perdomo  
**Institución:** Universidad Internacional de La Rioja (UNIR)  
**Facultad:** Facultad de Ciencias de la Salud  
**Programa:** Máster Universitario en Bioinformática  

