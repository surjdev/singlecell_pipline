#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
   Smart-seq2 Single-Cell RNA-seq Upstream Pipeline (Production DSL2)
========================================================================================
*/

include { FASTQC as FASTQC_RAW   } from './modules/fastqc'
include { FASTP                  } from './modules/fastp'
include { FASTQC as FASTQC_CLEAN } from './modules/fastqc'
include { STAR_ALIGN             } from './modules/star'
include { FEATURECOUNTS          } from './modules/featurecounts'
include { MULTIQC                } from './modules/multiqc'

workflow {
    log.info """
    ======================================================================
    🧬 PRODUCTION SMART-SEQ2 SINGLE-CELL UPSTREAM NEXTFLOW PIPELINE
    ======================================================================
    Samplesheet   : ${params.input ?: 'None (Using FASTQ glob)'}
    Reads pattern : ${params.reads ?: 'None'}
    STAR index    : ${params.star_index ?: 'NOT SPECIFIED'}
    GTF annotation: ${params.gtf ?: 'NOT SPECIFIED'}
    Output dir    : ${params.outdir}
    Min length    : ${params.min_length} bp
    Strandedness  : ${params.strandedness}
    ======================================================================
    """

    // 1. Validate required reference inputs
    if (!params.star_index) {
        error "ERROR: Please specify --star_index /path/to/star_index (or configure nextflow.config)"
    }
    if (!params.gtf) {
        error "ERROR: Please specify --gtf /path/to/genes.gtf (or configure nextflow.config)"
    }

    ch_star_index = file(params.star_index, checkIfExists: true)
    ch_gtf        = file(params.gtf, checkIfExists: true)

    // 2. Channel initialization: from samplesheet CSV or FASTQ filepairs
    if (params.input) {
        Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def meta = row.sample
                def r1 = file(row.fastq_1, checkIfExists: true)
                def r2 = file(row.fastq_2, checkIfExists: true)
                return tuple(meta, [r1, r2])
            }
            .set { ch_reads }
    } else if (params.reads) {
        Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .set { ch_reads }
    } else {
        error "ERROR: Please provide either --input samplesheet.csv or --reads 'data/raw/*_R{1,2}*.fastq.gz'"
    }

    // 3. Raw FastQC
    FASTQC_RAW(ch_reads)

    // 4. Quality & Adapter Trimming (fastp)
    FASTP(ch_reads)

    // 5. Post-Trimming FastQC
    FASTQC_CLEAN(FASTP.out.reads)

    // 6. Splice-Aware STAR Alignment & BAM Indexing
    STAR_ALIGN(FASTP.out.reads, ch_star_index)

    // 7. Subread featureCounts Gene Quantification
    ch_all_bams = STAR_ALIGN.out.bam.map { meta, bam, bai -> bam }.collect()
    FEATURECOUNTS(ch_all_bams, ch_gtf)

    // 8. MultiQC Summary Report
    ch_multiqc_files = Channel.empty()
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_CLEAN.out.zip)
        .mix(FASTP.out.json)
        .mix(STAR_ALIGN.out.log_final)
        .mix(STAR_ALIGN.out.gene_counts)
        .mix(FEATURECOUNTS.out.summary)
        .collect()

    MULTIQC(ch_multiqc_files)
}
