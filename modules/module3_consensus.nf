process GENERATE_CONSENSUS {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(filtered_fastq), path(assignments)

    output:
    tuple val(sample_id), path("${sample_id}_clusters_consensus.fasta"), emit: fasta

    script:
    """
    module3_consensus.py \\
        --sample $sample_id \\
        --fastq $filtered_fastq \\
        --assignments $assignments \\
        --racon_rounds 4
    """
}
