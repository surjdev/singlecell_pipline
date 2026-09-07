/*
========================================================================================
   Nextflow Module: MultiQC Report Aggregation
========================================================================================
*/

process MULTIQC {
    tag "multiqc"
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("single_cell_multiqc_report.html"), emit: report
    path("single_cell_multiqc_report_data"), emit: data

    script:
    """
    multiqc \\
        --filename "single_cell_multiqc_report.html" \\
        --force \\
        --interactive \\
        .
    """
}
