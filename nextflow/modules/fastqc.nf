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
    def stage = task.process.endsWith('FASTQC_CLEAN') ? 'clean' : 'raw'
    """
    ln -s '${reads[0]}' '${meta_id}_${stage}_R1.fastq.gz'
    ln -s '${reads[1]}' '${meta_id}_${stage}_R2.fastq.gz'
    fastqc \
        --threads ${task.cpus} \
        --noextract \
        --quiet \
        '${meta_id}_${stage}_R1.fastq.gz' '${meta_id}_${stage}_R2.fastq.gz'
    """
}
