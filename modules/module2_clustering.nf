process CLUSTERING {
    tag "$sample_id"
    label 'process_high'

    input:
    // Recibe el ID y el path del FASTQ filtrado del Módulo 1
    tuple val(sample_id), path(filtered_fastq)

    output:
    // 1. Asignaciones de cada lectura a un cluster
    tuple val(sample_id), path("${sample_id}_read_cluster_assignments.tsv"), emit: assignments
    
    // 2. Tabla de abundancias (OTU table)
    path "${sample_id}_OTU_table.tsv", emit: otu_table
    
    // 3. Gráfico de eficiencia (Ruido vs Clustered)
    path "${sample_id}_clustering_efficiency_plot.png", emit: plot
    
    // 4. Resumen numérico de los clusters
    path "${sample_id}_cluster_summary.tsv", emit: summary

    script:
    // Invocamos el script de Python con los parámetros definidos en params.config
    """
    module2_clustering.py \\
        --input ${filtered_fastq} \\
        --sample ${sample_id} \\
        --min_len ${params.min_length} \\
        --max_len ${params.max_length} \\
        --pca_comp ${params.pca_components}
    """
}
