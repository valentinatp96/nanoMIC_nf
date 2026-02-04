process VISUALIZATIONS {
    label 'process_medium'
    publishDir "${params.outdir}/module5_visualizations", mode: 'copy'

    input:
    path otu_tax_table // Recibe el archivo 'OTU_table_taxonomy.tsv'

    output:
    path "plots/*.png", emit: plots

    script:
    """
    # Creamos la carpeta de salida
    mkdir -p plots
    
    # Ejecutamos el script de python pasando el archivo de entrada y la carpeta de salida
    python3 ${baseDir}/bin/module5_visualizations.py ${otu_tax_table} plots
    """
}
