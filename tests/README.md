# Repeatable verification

From the repository root, with Pixi installed:

```bash
python tests/make_smoke_data.py /tmp/sc-demo
pixi run STAR --runMode genomeGenerate \
  --genomeDir /tmp/sc-demo/index --genomeFastaFiles /tmp/sc-demo/genome.fa \
  --sjdbGTFfile /tmp/sc-demo/genes.gtf --sjdbOverhang 74 \
  --genomeSAindexNbases 5 --runThreadN 2 --outFileNamePrefix /tmp/sc-demo/index/
pixi run bash scripts/run_all.sh --indir /tmp/sc-demo \
  --star-index /tmp/sc-demo/index --gtf /tmp/sc-demo/genes.gtf \
  --outdir /tmp/sc-demo/bash-results --threads 2
pixi run nextflow run nextflow/main.nf \
  --input /tmp/sc-demo/samples.csv --star_index /tmp/sc-demo/index \
  --gtf /tmp/sc-demo/genes.gtf --outdir /tmp/sc-demo/nf-results \
  --threads 2 -profile local -c nextflow/smoke.config \
  -work-dir /tmp/sc-demo/work
python tests/verify_smoke.py /tmp/sc-demo
# Assert identical count matrix from both execution routes:
cmp /tmp/sc-demo/bash-results/data/counts/gene_cell_count_matrix.tsv \
     /tmp/sc-demo/nf-results/counts/gene_cell_count_matrix.tsv
```

Expected: gene1 has 200 fragments each for cell_a and cell_b. Both routes must produce BAM/BAI,
FastQC raw/clean reports, fastp reports, featureCounts summary and MultiQC HTML/data.
This small, single-exon synthetic reference tests wiring, not splice junction sensitivity or real-data accuracy.
Run Nextflow again with `-resume` to check caching. Use a new output directory when changing inputs in Bash.

Python regressions:

```bash
cd workbench
pixi run test
```

Open `workbench/notebooks/00_end_to_end.ipynb` in the workbench Python kernel and run all cells.
The notebook uses synthetic expression by default and writes to `workbench/results/tutorial/`.
