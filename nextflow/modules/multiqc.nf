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
    path(qc_files), stageAs: 'inputs_??????/*'

    output:
    path("single_cell_multiqc_report.html"), emit: report
    path("single_cell_multiqc_report_data"), emit: data

    script:
    def interactive_flag = params.interactive_plots ? "--interactive" : ""
    """
    multiqc \\
        --filename "single_cell_multiqc_report.html" \\
        --force \\
        --fullnames \\
        ${interactive_flag} \\
        .
    """
}
