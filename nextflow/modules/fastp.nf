/*
========================================================================================
   Nextflow Module: fastp Quality & Adapter Trimming (Production Ready)
========================================================================================
*/

process FASTP {
    tag "${meta_id}"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta_id), path(reads)

    output:
    tuple val(meta_id), path("${meta_id}_R{1,2}.clean.fastq.gz"), emit: reads
    path("${meta_id}_fastp.json")                                , emit: json
    path("${meta_id}_fastp.html")                                , emit: html

    script:
    def adapter_args = "--detect_adapter_for_pe"
    if (params.adapter_r1 != 'auto' && params.adapter_r1 != '') {
        adapter_args = "--adapter_sequence ${params.adapter_r1}"
        if (params.adapter_r2 != 'auto' && params.adapter_r2 != '') {
            adapter_args += " --adapter_sequence_r2 ${params.adapter_r2}"
        }
    }
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta_id}_R1.clean.fastq.gz \\
        --out2 ${meta_id}_R2.clean.fastq.gz \\
        ${adapter_args} \\
        --trim_poly_g \\
        --poly_g_min_len 10 \\
        --trim_poly_x \\
        --poly_x_min_len 10 \\
        --cut_front \\
        --cut_front_window_size 4 \\
        --cut_front_mean_quality ${params.min_quality} \\
        --cut_tail \\
        --cut_tail_window_size 4 \\
        --cut_tail_mean_quality ${params.min_quality} \\
        --qualified_quality_phred ${params.min_quality} \\
        --length_required ${params.min_length} \\
        --thread ${task.cpus} \\
        --json ${meta_id}_fastp.json \\
        --html ${meta_id}_fastp.html
    """
}
