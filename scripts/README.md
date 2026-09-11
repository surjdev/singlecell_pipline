# คู่มือใช้ scripts เพื่อเรียนรู้และประมวลผลข้อมูล

ทุกคำสั่งรันจาก repository root ภายใน `pixi run ...` หรือ `pixi shell`
script แต่ละตัวเป็น Bash wrapper อ่านตัวเลือกจริงของเครื่องมือได้ใน source
ดูคู่มือการรันฉบับเต็มได้ที่ [คู่มือการรันฉบับละเอียด (RUN_GUIDE)](../docs/RUN_GUIDE.md) และดูความหมายเชิงละเอียดของ flags ที่ [TOOLS_RUNBOOK](../docs/TOOLS_RUNBOOK.md)

| Script | Input → output | สิ่งที่ควรตรวจ |
|---|---|---|
| `setup_reference.sh` | FASTA + GTF → STAR index | species/build ตรงกัน; index ใช้ STAR รุ่นที่เข้ากัน |
| `qc_data_00.sh` | FASTQ → HTML/ZIP | base quality, adapters, read length; warning ไม่ใช่เหตุให้ทิ้งเซลล์ทันที |
| `fastp_01.sh` | gzip paired FASTQ → clean FASTQ + JSON/HTML | reads retained, adapter trimming, length distribution |
| `align_02.sh` | clean FASTQ + index → BAM/BAI + STAR logs | uniquely mapped %, multimapping, unmapped reasons |
| `quant_03.sh` | BAM + GTF → counts TSV + assignment summary | Assigned %, gene IDs, strandedness |
| `multiqc_04.sh` | QC/log directories → HTML dashboard | เปรียบเทียบเซลล์และ raw/clean reports |
| `run_all.sh` | raw FASTQ + references → ผลครบทุก stage | ใช้ output directory ใหม่ต่อชุดข้อมูล |

## รันทีละขั้นตอน

```bash
pixi run bash scripts/qc_data_00.sh data/raw reports/qc_raw 4
pixi run bash scripts/fastp_01.sh data/raw data/clean reports/fastp 4
pixi run bash scripts/qc_data_00.sh data/clean reports/qc_clean 4
pixi run bash scripts/align_02.sh data/clean data/aligned data/reference/star_index 4
pixi run bash scripts/quant_03.sh data/aligned data/counts data/reference/genes.gtf 4
pixi run bash scripts/multiqc_04.sh reports/multiqc reports/qc_raw reports/qc_clean reports/fastp data/aligned data/counts
```

positional arguments เรียงตามคำสั่งข้างต้น; ถ้าไม่ใส่จะใช้ค่าใน `scripts/env.sh`
`run_all.sh --help` และ `setup_reference.sh --help` แสดง named options

## การตั้งค่า

คัดลอก `config.env.example` เป็น `config.env` แล้วปรับค่าตามงาน หรือส่ง `--config FILE` ให้ runner
ไฟล์ config เป็น shell code ควรใช้ไฟล์ที่เชื่อถือได้
runner โหลด config ก่อน แล้วใช้ CLI override; ค่าที่ resolve แล้วส่งให้ child stages
`--outdir DIR` สร้าง `DIR/data/{clean,aligned,counts}` และ `DIR/reports/*`
Bash ไม่มี cache/resume อัตโนมัติ: การรันซ้ำเขียนทับไฟล์ และอาจรวม BAM เก่าที่ค้างอยู่ใน directory เดิม
ใช้ปลายทางใหม่ต่อ dataset หรือใช้ Nextflow `-resume`

| ค่า config | ค่าเริ่มต้น | ผลต่อการประมวลผล |
|---|---:|---|
| THREADS | CPU ที่ตรวจพบ | ควรตั้ง 4–8 เอง; fastp ถูกจำกัดไม่เกิน 16 |
| MIN_READ_LENGTH | 35 | สั้นกว่านี้หลัง trim จะถูกทิ้ง |
| MIN_QUALITY | 20 | quality threshold ของ fastp |
| MAX_UNQUALIFIED_PCT | 30 | สัดส่วน base คุณภาพต่ำที่ยอมรับ |
| ADAPTER_FWD / ADAPTER_REV | auto | ระบุ adapter ตาม library prep หากทราบ |
| STRANDEDNESS | 0 | featureCounts: 0 unstranded, 1 forward, 2 reverse |
| FEATURE_TYPE / ATTRIBUTE_TYPE | exon / gene_id | รวม reads ตาม exon เป็น gene counts |
| MIN_MAPQ | 10 | mapping quality ต่ำกว่านี้ไม่ถูกนับ |
| STAR_RAM_LIMIT | 31000000000 | sort buffer bytes ไม่ใช่เพดาน RAM ของ process |

fastp ปัจจุบันเปิด poly-G/poly-X และ front/tail trimming จึงควรตรวจ retained reads กับ library ของตน
ไม่ได้ทำ duplicate removal: reads ซ้ำอาจเกิดจาก expression จริง โดยเฉพาะ protocol ที่ไม่มี UMI

## Naming และ reference

รองรับ `cell_R1.fastq.gz/cell_R2.fastq.gz`, `_R1_001/_R2_001`, `.R1/.R2`, `_1/_2`, และ `_R1.fq.gz/_R2.fq.gz`
ค้นหาไม่เกิน 2 ระดับ directory และรองรับ symlink; ห้ามมี sample basename ซ้ำข้าม directory
lane แยกถูกมองเป็นคนละ sample: รวม lanes ต่อ cell ก่อนรัน
FastQC รับ uncompressed FASTQ ได้ แต่ fastp/STAR wrappers นี้ใช้ gzip

```bash
pixi run bash scripts/setup_reference.sh \
  --fasta /refs/genome.fa --gtf /refs/genes.gtf \
  --outdir data/reference --threads 8 --sjdbOverhang 74
```

custom FASTA/GTF ถูกใช้งานจาก path เดิม ไม่ได้คัดลอกเป็น `data/reference/genes.gtf`
ตอน quantify จึงต้องส่ง `--gtf /refs/genes.gtf` เดิม
`sjdbOverhang` ใช้ read length − 1; genome เล็กสำหรับการทดสอบต้องปรับ `genomeSAindexNbases` ด้วยคำสั่ง STAR โดยตรงตาม [smoke test](../tests/README.md)

## Output contract สำหรับ downstream

`gene_cell_count_matrix.tsv` มี header `Geneid` แล้วชื่อเซลล์ แถวคือ gene ID
`featurecounts_raw.txt` เก็บรายละเอียด annotation และ `.summary` เก็บสาเหตุ reads ที่ไม่ถูกนับ
ใช้ GTF รุ่นเดียวกันเติม gene symbols ให้ `adata.var['gene_name']` ถ้า IDs เป็น Ensembl;
การตรวจ prefix `MT-` บน Ensembl IDs อย่างเดียวจะไม่พบ mitochondrial genes
