process QC_READS {
    tag "$sample_id"
    label 'env_nanopore' // Etiqueta opcional para organización
    
    // Definimos etiquetas para el manejo de recursos en el config
    label 'process_medium'

    input:
    tuple val(sample_id), path(sample_dir)

    output:
    // 1. Las lecturas filtradas para el siguiente módulo (Clustering)
    tuple val(sample_id), path("${sample_id}_V4_filtered.fastq.gz"), emit: reads
    
    // 2. Archivos para visualización y MultiQC
    path "nanoplot_${sample_id}", emit: nanoplot_logs
    path "${sample_id}_cutadapt.log", emit: cutadapt_logs
    path "${sample_id}_concatenated.fastq.gz", emit: concat_fastq

    script:
    // Usamos el script de bin/ aprovechando que NF lo pone en el PATH
    """
    module1_qc.sh \\
        --sample "${sample_id}" \\
        --input_dir "${sample_dir}" \\
        --primer_f "${params.primer_f}" \\
        --primer_r "${params.primer_r}"
    """
}
