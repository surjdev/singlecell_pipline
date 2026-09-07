/*
========================================================================================
   Nextflow Module: FastQC Quality Control
========================================================================================
*/

process FASTQC {
    tag "${meta_id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: '*.{html,zip}'

    input:
    tuple val(meta_id), path(reads)

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --noextract \\
        --quiet \\
        ${reads}
    """
}
