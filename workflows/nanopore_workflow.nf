nextflow.enable.dsl = 2

include { MODULE1_QC } from '../modules/module1_qc'

workflow NANOPORE_WORKFLOW {

    take:
    samples

    main:
    qc = MODULE1_QC(samples)

    emit:
    qc.filtered_reads
}
