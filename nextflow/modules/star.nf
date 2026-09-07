/*
========================================================================================
   Nextflow Module: STAR Splice-Aware Alignment (Production Ready)
========================================================================================
*/

process STAR_ALIGN {
    tag "${meta_id}"
    label 'process_high'
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: '*.{bam,bai,out,tab}'

    input:
    tuple val(meta_id), path(reads)
    path(index)

    output:
    tuple val(meta_id), path("${meta_id}_Aligned.sortedByCoord.out.bam"), path("${meta_id}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    path("${meta_id}_Log.final.out")                                                                                           , emit: log_final
    path("${meta_id}_ReadsPerGene.out.tab")                                                                                    , emit: gene_counts

    script:
    def sort_ram = task.memory ? (task.memory.toBytes() * 0.75).longValue() : 31000000000
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMattributes NH HI AS nM NM MD jM jI XS \\
        --outFilterType BySJout \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 10 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --limitBAMsortRAM ${sort_ram} \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${meta_id}_

    samtools index -@ ${task.cpus} ${meta_id}_Aligned.sortedByCoord.out.bam
    """
}
