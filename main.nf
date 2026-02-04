#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/* ============================
    1. Importar Módulos
   ============================ */
include { QC_READS } from './modules/module1_qc'
include { CLUSTERING } from './modules/module2_clustering'
include { GENERATE_CONSENSUS } from './modules/module3_consensus'
include { TAXONOMY } from './modules/module4_taxonomy'
include { VISUALIZATIONS } from './modules/module5_visualizations'
include { MULTIQC } from './modules/module6_multiqc'


/* ============================
    2. Workflow Principal
   ============================ */
workflow {
    // 1. Crear canal de entrada
    ch_samples = Channel
        .fromPath("${params.input_dir}/*", type: 'dir', checkIfExists: true)
        .map { dir -> tuple(dir.name, dir) }

    // 2. Ejecutar Módulo 1: QC
    QC_READS(ch_samples)
    
    // 3. Ejecutar Módulo 2: Clustering
    CLUSTERING(QC_READS.out.reads)

    // 4. Ejecutar Módulo 3: Consenso
    // Combinamos el FASTQ filtrado con las asignaciones de clusters
    ch_for_consensus = QC_READS.out.reads.join(CLUSTERING.out.assignments)
    GENERATE_CONSENSUS(ch_for_consensus)

    // 5. Ejecutar Módulo 4: Taxonomía (Global)
    TAXONOMY(
    GENERATE_CONSENSUS.out.fasta.map{ it[1] }.collect(),
        CLUSTERING.out.assignments.map{ it[1] }.collect(),
        file(params.silva_seqs),
        file(params.silva_tax)
    )
    
    // 6. Ejecutar Módulo 5: Visualizaciones
    // Este toma la tabla de taxonomía final generada en el paso anterior
    VISUALIZATIONS(TAXONOMY.out.otu_tax_table)

// 7. Módulo 6: MultiQC (Global)
    // Usamos los nombres reales que definiste en el módulo QC_READS
    ch_multiqc_files = QC_READS.out.nanoplot_logs.collect()
        .mix(QC_READS.out.cutadapt_logs.collect())
    
    MULTIQC(ch_multiqc_files)
}

/* ============================
    Mensaje de finalización
   ============================ */
workflow.onComplete {
    println """
    🚀 Pipeline nanoMIC finalizado!
    🏁 Estado: ${workflow.success ? 'OK' : 'FALLÓ'}
    📂 Resultados en: ${params.outdir}
    """
}
