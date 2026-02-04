process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/module1_qc/multiqc_report", mode: 'copy'

    input:
    path ('logs/*') // Recibe una lista de archivos de log/stats

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc logs/ --filename multiqc_report.html
    """
}
