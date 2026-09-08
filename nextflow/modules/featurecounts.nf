/*
========================================================================================
   Nextflow Module: Subread featureCounts Quantification (Production Ready)
========================================================================================
*/

process FEATURECOUNTS {
    tag "all_samples"
    label 'process_medium'
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: '*.{tsv,txt,summary}'

    input:
    path(bams)
    path(gtf)

    output:
    path("gene_cell_count_matrix.tsv")   , emit: matrix
    path("featurecounts_raw.txt")        , emit: raw_counts
    path("featurecounts_raw.txt.summary"), emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p \\
        --countReadPairs \\
        -t exon \\
        -g gene_id \\
        -a '${gtf}' \\
        -o featurecounts_raw.txt \\
        -s ${params.strandedness} \\
        -Q 10 \\
        ${bams.collect { "\"${it}\"" }.join(" ")}

    # Format clean Gene x Cell expression matrix
    awk '
    BEGIN { FS="\\t"; OFS="\\t" }
    /^#/ { next }
    \$1 == "Geneid" {
        printf "%s", \$1
        for (i=7; i<=NF; i++) {
            n = split(\$i, parts, "/")
            fname = parts[n]
            sub(/_Aligned\\.sortedByCoord\\.out\\.bam\$/, "", fname)
            printf "\\t%s", fname
        }
        printf "\\n"
        next
    }
    {
        printf "%s", \$1
        for (i=7; i<=NF; i++) {
            printf "\\t%s", \$i
        }
        printf "\\n"
    }
    ' featurecounts_raw.txt > gene_cell_count_matrix.tsv
    """
}
