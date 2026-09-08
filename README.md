# Smart-seq2 processing and downstream analysis

สำหรับ paired-end Smart-seq2 / plate-based RNA-seq: หนึ่งคู่ FASTQ = หนึ่งเซลล์
Bash และ Nextflow ทำ FastQC → fastp → FastQC → STAR/samtools → featureCounts → MultiQC
ผลลัพธ์เป็น **fragment counts** (genes × cells) ไม่ใช่ UMI matrix และไม่ใช่ workflow สำหรับ raw 10x FASTQ

## เริ่มใช้งาน

```bash
pixi install
pixi run bash scripts/setup_reference.sh --species human --threads 8
pixi run bash scripts/run_all.sh \
  --indir data/raw --star-index data/reference/star_index \
  --gtf data/reference/genes.gtf --outdir results/bash --threads 4
```

Human reference/index ต้องใช้ RAM และพื้นที่มาก โปรดกำหนดทรัพยากรตาม reference;
`STAR_RAM_LIMIT` เป็น sort buffer เพิ่มจากหน่วยความจำ genome ไม่ใช่ RAM รวม
ศึกษาวิธีรันแต่ละ stage และ parameter ใน [คู่มือ scripts](scripts/README.md)
และ [รายละเอียดเครื่องมือ](docs/TOOLS_RUNBOOK.md)

## Nextflow

ตัวอย่าง workflow ที่รันได้อยู่ใน [nextflow/main.nf](nextflow/main.nf)

```bash
pixi run nextflow run nextflow/main.nf \
  --input nextflow/samplesheet.example.csv \
  --star_index /absolute/path/star_index --gtf /absolute/path/genes.gtf \
  --outdir results/nextflow --threads 4 -profile local
# เพิ่ม -resume เพื่อใช้ task ที่สำเร็จแล้วจาก cache
```

แก้ samplesheet ก่อนรัน: `sample,fastq_1,fastq_2`; sample ต้องไม่ซ้ำ
relative FASTQ paths อ้างจาก directory ของ CSV ส่วน reference paths อ้างจาก working directory
ใช้ `.fastq.gz` / `.fq.gz` และชื่อ sample เป็นตัวอักษร ตัวเลข `_ . -`
อ่าน [คู่มือ Nextflow](nextflow/README.md) สำหรับทรัพยากรและข้อจำกัด profiles

## Downstream บน Jupyter

```bash
cd workbench
pixi run lab
```

เริ่มที่ [00_end_to_end.ipynb](workbench/notebooks/00_end_to_end.ipynb)
รันครบได้ด้วยข้อมูลสังเคราะห์ ไม่ดาวน์โหลดข้อมูล แล้วเปลี่ยน `COUNTS_PATH` เพื่อรับผล upstream:

- Bash: `results/bash/data/counts/gene_cell_count_matrix.tsv`
- Nextflow: `results/nextflow/counts/gene_cell_count_matrix.tsv`

[sc_workbench](workbench/README.md) ใช้ AnnData เป็นข้อมูลหลัก
ใช้ Pandas กับ `adata.obs/var`, NumPy/SciPy กับ `adata.X/layers`, Scikit-learn กับ `adata.obsm['X_pca']`
fluent API เป็นทางเลือก ไม่จำเป็นต้องใช้ wrapper สำหรับทุกขั้นตอน

## ตรวจสอบซ้ำ

ดู [ผลตรวจคุณภาพและขอบเขตการทดสอบ](docs/QUALITY_REVIEW.md) และ
[วิธี smoke test ด้วยเครื่องมือจริง](tests/README.md)
การทดสอบสังเคราะห์ยืนยันการต่อระบบ ไม่ได้ยืนยันคุณภาพ mapping/annotation บนข้อมูลทดลองจริง
