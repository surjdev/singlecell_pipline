# Nextflow execution guide

`main.nf` is the runnable DSL2 example, sharing no execution dependency with the Bash wrappers.
Run from the repository root using `pixi run nextflow run nextflow/main.nf ...`.

```bash
pixi run nextflow run nextflow/main.nf --help
pixi run nextflow run nextflow/main.nf \
  --input samples.csv --star_index /refs/star_index --gtf /refs/genes.gtf \
  --outdir results --threads 4 -profile local -resume
```

CSV columns are `sample,fastq_1,fastq_2`. Relative read paths resolve against the CSV directory.
All samples are validated before tasks launch: empty/missing fields, duplicate IDs, reused read paths,
missing files and identical mates fail. One gzip paired-end library corresponds to one cell.
Use simple file names without quotes or shell metacharacters. Merge lanes before preparing the sheet.
The alternative `--reads 'data/raw/*_R{1,2}.fastq.gz'` uses Nextflow's filename pairing.

Output: `fastqc/`, `fastp/`, `aligned/`, `counts/`, `multiqc/`, `pipeline_info/`.
Clean FASTQ stays in the work directory; retain it for `-resume`. Count columns are sorted by BAM name.
MultiQC stages reports in separate directories to avoid input basename collisions.

`--threads` caps CPUs per process (low: 2, medium: 6, high: 12).
STAR's default process allocation is 36 GB; `--star_sort_ram` is a separate sort buffer (4 GB default).
The genome index itself needs additional RAM. Concurrent tasks increase total machine memory use.
Override resources with a config passed by `-c`:

```groovy
process {
    withName: STAR_ALIGN {
        cpus = 8
        memory = 48.GB
        maxForks = 1
    }
}
params.star_sort_ram = 4000000000
```

`smoke.config` is exclusively for the tiny synthetic reference, not human/mouse runs.
The `local` profile with installed Pixi tools has been exercised. Docker/Singularity image availability
and Slurm execution have not been tested here; verify the configured images and cluster queue before use.
The former empty `conda` profile was removed because it specified no tool environments.
Unused `max_memory`, `max_cpus`, `max_time` placeholders were removed;
use process config overrides as above. Retries currently repeat resource allocations without increasing RAM.

See [comprehensive execution guide](../docs/RUN_GUIDE.md) and [repeatable real-tool smoke test](../tests/README.md).
Nextflow staging semantics: [official process documentation](https://nextflow.io/docs/stable/process.html).
