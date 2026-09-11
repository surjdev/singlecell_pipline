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
        log.info """Usage: nextflow run nextflow/main.nf --sample_sheet samples.csv --star_index INDEX --gtf genes.gtf --outdir results -profile local
CSV columns: sample,fastq_1[,fastq_2]. Single-End (fastq_1 only) or Paired-End (fastq_1 and fastq_2) are supported.
Relative FASTQ paths resolve against the CSV directory.
One unique sample per cell; gzip reads required. Use -resume to reuse successful tasks."""
        return
    }
    if (!(params.strandedness.toString() in ['0','1','2'])) error 'strandedness must be 0, 1 or 2'

    def sample_sheet = params.sample_sheet ?: params.samplesheet ?: params.input

    log.info """
    ======================================================================
    🧬 PRODUCTION SMART-SEQ2 SINGLE-CELL UPSTREAM NEXTFLOW PIPELINE
    ======================================================================
    Samplesheet   : ${sample_sheet ?: 'None (Using FASTQ glob)'}
    Reads pattern : ${params.reads ?: 'None'}
    Single-End    : ${params.single_end != null ? params.single_end : 'Auto-detect'}
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
    // 2. Channel initialization: from samplesheet CSV or FASTQ filepairs/single-end
    if (sample_sheet) {
        def sheetDir = file(sample_sheet).toAbsolutePath().parent
        Channel
            .fromPath(sample_sheet, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def meta = row.sample ?: row.sample_id
                def fq1  = row.fastq_1 ?: row.fq1 ?: row.read1 ?: row.fastq
                def fq2  = row.fastq_2 ?: row.fq2 ?: row.read2 ?: null

                if (!meta || !fq1) error "CSV requires at least 'sample' and 'fastq_1' columns"
                if (!(meta ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample ID: ${meta}"

                def r1 = file(sheetDir.resolve(fq1).normalize(), checkIfExists: true)
                if (!r1.name.endsWith('.gz')) error "Expected gzip compressed read: ${r1}"

                def reads = [r1]
                if (fq2 && fq2.trim() != '') {
                    def r2 = file(sheetDir.resolve(fq2).normalize(), checkIfExists: true)
                    if (!r2.name.endsWith('.gz')) error "Expected gzip compressed read: ${r2}"
                    if (r1 == r2) error "Expected distinct gzip mates for sample: ${meta}"
                    reads.add(r2)
                }
                return tuple(meta, reads)
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
        if (params.single_end) {
            Channel
                .fromPath(params.reads, checkIfExists: true)
                .map { f ->
                    def sample_id = f.name.replaceAll(/(\.fastq|\.fq)?\.gz$/, '')
                    return tuple(sample_id, [f])
                }
                .set { ch_reads }
        } else {
            Channel
                .fromFilePairs(params.reads, checkIfExists: true)
                .set { ch_reads }
        }
    } else {
        error "ERROR: Please provide either --sample_sheet samplesheet.csv or --reads 'data/raw/*_R{1,2}*.fastq.gz'"
    }

    // Determine whether single-end or paired-end for quantification
    ch_is_single_end = ch_reads.first().map { meta, reads ->
        params.single_end != null ? (params.single_end as boolean) : (reads.size() == 1)
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
    FEATURECOUNTS(ch_all_bams, ch_gtf, ch_is_single_end)

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
