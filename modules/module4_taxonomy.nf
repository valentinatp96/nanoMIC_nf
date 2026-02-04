process TAXONOMY {
    label 'process_high'
    // El entorno conda se activa automáticamente gracias a tu nextflow.config
    
    publishDir "${params.outdir}/module4_taxonomy", mode: 'copy'

    input:
    path fastas              // Lista de archivos .fasta de todas las muestras
    path assignments         // Lista de archivos .tsv de todas las muestras
    path db_seqs             // Archivo silva-xxx-seqs.qza
    path db_tax              // Archivo silva-xxx-tax.qza

    output:
    path "OTU_table_taxonomy.tsv", emit: otu_tax_table
    path "taxonomy_vsearch.tsv",    emit: taxonomy_tsv
    path "consensus_renamed.fasta", emit: renamed_fasta
    path "taxonomy_assignment_qc.png", emit: qc_plot
    path "*.qza",                   emit: qza_files

    script:
    """
    # 1. Consolidar inputs de múltiples muestras en archivos únicos
    cat ${fastas} > all_consensus_merged.fasta
    
    # Unir asignaciones de clusters omitiendo cabeceras repetidas
    awk 'FNR==1 && NR!=1{next;}{print}' ${assignments} > all_assignments_merged.tsv

    # 2. Ejecutar el script enviando los 4 argumentos requeridos
    bash module4_taxonomy.sh \\
        all_consensus_merged.fasta \\
        all_assignments_merged.tsv \\
        ${db_seqs} \\
        ${db_tax}
    """
}
