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
    if (params.help) {
        log.info """Usage: nextflow run nextflow/main.nf --input samples.csv --star_index INDEX --gtf genes.gtf --outdir results -profile local
CSV columns: sample,fastq_1,fastq_2. Relative FASTQ paths resolve against the CSV directory.
One unique sample per cell; gzip paired reads required. Use -resume to reuse successful tasks."""
        return
    }
    if (!(params.strandedness.toString() in ['0','1','2'])) error 'strandedness must be 0, 1 or 2'

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

    ['Genome', 'SA', 'SAindex'].each { name ->
        if (!ch_star_index.resolve(name).exists()) error "Incomplete STAR index: missing ${name}"
    }
    // 2. Channel initialization: from samplesheet CSV or FASTQ filepairs
    if (params.input) {
        def sheetDir = file(params.input).toAbsolutePath().parent
        Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                if (!row.sample || !row.fastq_1 || !row.fastq_2) error 'CSV requires sample,fastq_1,fastq_2'
                def meta = row.sample
                if (!(meta ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample ID: ${meta}"
                def r1 = file(sheetDir.resolve(row.fastq_1).normalize(), checkIfExists: true)
                def r2 = file(sheetDir.resolve(row.fastq_2).normalize(), checkIfExists: true)
                if (r1 == r2 || !r1.name.endsWith('.gz') || !r2.name.endsWith('.gz')) error "Expected distinct gzip mates: ${meta}"
                return tuple(meta, [r1, r2])
            }
            .toList()
            .map { rows ->
                if (!rows) error 'Samplesheet is empty'
                if (rows.collect { it[0] }.unique().size() != rows.size()) error 'Duplicate sample IDs'
                def paths = rows.collectMany { it[1] }
                if (paths.unique(false).size() != paths.size()) error 'FASTQ files reused across samples'
                return rows
            }
            .flatMap { it }
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
    ch_all_bams = STAR_ALIGN.out.bam.map { meta, bam, bai -> bam }.collect(sort: true)
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
