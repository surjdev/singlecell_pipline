# 🚀 คู่มือการรัน Single-Cell RNA-seq Pipeline (Bash Scripts & Nextflow) ฉบับละเอียด

คู่มือนี้รวบรวมขั้นตอนการรันไปป์ไลน์วิเคราะห์ข้อมูล Smart-seq2 / Plate-based scRNA-seq แบบต้นน้ำ (Upstream) ตั้งแต่ไฟล์ดิบ FASTQ จนถึง **Gene-by-Cell Count Matrix** และ **MultiQC Report** อย่างละเอียด ทั้งแบบ **Pure-Bash Scripts** และ **Nextflow DSL2**

---

## 📑 สารบัญ
1. [ความต้องการของระบบและทรัพยากร (System Requirements)](#1-ความต้องการของระบบและทรัพยากร-system-requirements)
2. [การเตรียมสภาพแวดล้อมด้วย Pixi (Environment Setup)](#2-การเตรียมสภาพแวดล้อมด้วย-pixi-environment-setup)
3. [การเตรียม Reference Genome และสร้าง STAR Index](#3-การเตรียม-reference-genome-และสร้าง-star-index)
4. [Track 1: การรันด้วย Pure-Bash Scripts](#4-track-1-การรันด้วย-pure-bash-scripts)
   - [การตั้งค่า config.env](#41-การตั้งค่า-configenv)
   - [วิธีที่ 1: รันคำสั่งเดียวจบด้วย Master Runner (`run_all.sh`)](#42-วิธีที่-1-รันคำสั่งเดียวจบด้วย-master-runner-run_allsh)
   - [วิธีที่ 2: รันแยกทีละขั้นตอน (Step-by-Step)](#43-วิธีที่-2-รันแยกทีละขั้นตอน-step-by-step)
5. [Track 2: การรันด้วย Nextflow DSL2 Pipeline](#5-track-2-การรันด้วย-nextflow-dsl2-pipeline)
   - [การเตรียม Samplesheet (รองรับทั้ง Single-End และ Paired-End)](#51-การเตรียม-samplesheet)
   - [คำสั่งการรัน Nextflow พร้อม Options ที่สำคัญ](#52-คำสั่งการรัน-nextflow)
   - [การใช้ `-resume` เพื่อประหยัดเวลา](#53-การใช้--resume)
   - [Execution Profiles (Local, Docker, Singularity, Slurm)](#54-execution-profiles)
   - [การปรับแต่ง CPU/RAM แบบละเอียด](#55-การปรับแต่ง-cpuram-ใน-nextflow)
6. [โครงสร้างผลลัพธ์และการตรวจสอบคุณภาพ (Outputs & QC Interpretation)](#6-โครงสร้างผลลัพธ์และการตรวจสอบคุณภาพ)
7. [การส่งต่อไปยัง Downstream Analysis (Jupyter / AnnData)](#7-การส่งต่อไปยัง-downstream-analysis)
8. [ปัญหาที่พบบ่อยและวิธีแก้ไข (Troubleshooting & FAQs)](#8-ปัญหาที่พบบ่อยและวิธีแก้ไข)

---

## 1. ความต้องการของระบบและทรัพยากร (System Requirements)

- **ระบบปฏิบัติการ**: Linux (Ubuntu 20.04+, Debian, CentOS/RHEL)
- **RAM**:
  - ขั้นต่ำ **32 GB RAM** สำหรับ Human (GRCh38) หรือ Mouse (GRCm39) เนื่องจาก STAR index โหลดจีโนมเข้าหน่วยความจำ
  - หากใช้ Small Synthetic Genome (เช่น smoke test) ใช้เพียง 4–8 GB
- **Storage**:
  - พื้นที่ว่างอย่างน้อย 100–150 GB (ขึ้นอยู่กับจำนวนเซลล์ ขนาด STAR index ~30 GB + BAM files)
- **CPU**: แนะนำ 8–16 cores ขึ้นไป

---

## 2. การเตรียมสภาพแวดล้อมด้วย Pixi (Environment Setup)

โปรเจกต์นี้จัดการ bioinformatics tools ทั้งหมด (FastQC, fastp, STAR, samtools, Subread featureCounts, MultiQC, Nextflow) ผ่าน [Pixi](https://pixi.sh/):

### 2.1 ติดตั้ง Dependencies
รันจาก Root directory ของโปรเจกต์:
```bash
# ตรวจสอบการติดตั้ง dependencies
pixi install
```

### 2.2 การเรียกใช้งาน
มี 2 วิธีหลักในการรันคำสั่ง:
1. **ผ่านคำนำหน้า `pixi run`** (ไม่ต้อง activate shell):
   ```bash
   pixi run bash scripts/run_all.sh --help
   pixi run nextflow run nextflow/main.nf --help
   ```
2. **เข้าสู่ Environment Shell โดยตรง**:
   ```bash
   pixi shell
   # เมื่ออยู่ใน shell สามารถพิมพ์คำสั่งตรงๆ ได้ทันที:
   bash scripts/run_all.sh --help
   nextflow run nextflow/main.nf --help
   exit # เพื่อออกจาก environment
   ```

---

## 3. การเตรียม Reference Genome และสร้าง STAR Index

ก่อนรันข้อมูลจริง จำเป็นต้องมี **STAR Genome Index** และ **GTF Annotation File**

### กรณี A: ดาวน์โหลดและสร้าง Index อัตโนมัติ (Human / Mouse)
ใช้สคริปต์ `scripts/setup_reference.sh` ซึ่งจะดาวน์โหลดไฟล์จาก GENCODE / Ensembl และสร้าง STAR index ให้พร้อมใช้งาน:

```bash
# สำหรับมนุษย์ (Human GRCh38) - แนะนำใช้ RAM >= 32GB
pixi run bash scripts/setup_reference.sh \
  --species human \
  --outdir data/reference \
  --threads 16 \
  --sjdbOverhang 99

# สำหรับหนู (Mouse GRCm39)
pixi run bash scripts/setup_reference.sh \
  --species mouse \
  --outdir data/reference \
  --threads 16 \
  --sjdbOverhang 99
```
> [!NOTE]
> **สูตรการตั้งค่า `--sjdbOverhang`**: แนะนำให้ตั้งเป็น `Read Length - 1` เช่น
> - ความยาว Read 100 bp $\rightarrow$ `sjdbOverhang 99` (ค่าเริ่มต้นคือ 100)
> - ความยาว Read 75 bp $\rightarrow$ `sjdbOverhang 74`
> - ความยาว Read 150 bp $\rightarrow$ `sjdbOverhang 149`

### กรณี B: มีไฟล์ FASTA และ GTF อยู่แล้ว (Custom Reference)
หากมีไฟล์ `.fa` และ `.gtf` ของตนเองอยู่แล้ว สามารถสั่งสร้าง STAR Index ได้โดยไม่ต้องดาวน์โหลดใหม่:

```bash
pixi run bash scripts/setup_reference.sh \
  --fasta /path/to/my_genome.fa \
  --gtf /path/to/my_genes.gtf \
  --outdir data/reference \
  --threads 16 \
  --sjdbOverhang 99
```

ไฟล์ที่ได้ในโฟลเดอร์ปลายทาง (`data/reference/`):
- `genome.fa`: Genome FASTA
- `genes.gtf`: Gene annotation GTF
- `star_index/`: โฟลเดอร์ STAR binary index (`Genome`, `SA`, `SAindex`)

---

## 4. Track 1: การรันด้วย Pure-Bash Scripts

ระบบ Bash เป็น modular scripts แยกตามขั้นตอน เหมาะสำหรับการเรียนรู้ การทดสอบทีละสเต็ป หรือกรณีที่เครื่องไม่มี Nextflow

### 4.1 การตั้งค่า `config.env`

คัดลอกไฟล์ตัวอย่างเพื่อกำหนดพารามิเตอร์ส่วนตัว:
```bash
cp config.env.example config.env
```
เปิดแก้ไขไฟล์ `config.env` ตามต้องการ:
```bash
# กำหนด Path ข้อมูลและ Reference
RAW_DIR="data/raw"
CLEAN_DIR="data/clean"
ALIGNED_DIR="data/aligned"
COUNTS_DIR="data/counts"
REPORTS_DIR="reports"
STAR_INDEX="data/reference/star_index"
GTF_FILE="data/reference/genes.gtf"

# Hardware Resources
THREADS=8
STAR_RAM_LIMIT=31000000000 # Sort buffer (~31 GB)

# Filtering & Trimming (fastp)
MIN_READ_LENGTH=35
MIN_QUALITY=20
MAX_UNQUALIFIED_PCT=30
ADAPTER_FWD="CTGTCTCTTATACACATCT" # Nextera transposase adapter สำหรับ Smart-seq2 (หรือ "auto")
ADAPTER_REV="CTGTCTCTTATACACATCT"

# Quantification (featureCounts)
STRANDEDNESS=0 # 0 = unstranded (Smart-seq2 ส่วนใหญ่), 1 = stranded, 2 = reverse
FEATURE_TYPE="exon"
ATTRIBUTE_TYPE="gene_id"
MIN_MAPQ=10
```

---

### 4.2 วิธีที่ 1: รันคำสั่งเดียวจบด้วย Master Runner (`run_all.sh`)

สคริปต์ `scripts/run_all.sh` จะรันครบทั้ง 5 ขั้นตอน (FastQC ดิบ $\rightarrow$ fastp $\rightarrow$ FastQC คลีน $\rightarrow$ STAR $\rightarrow$ featureCounts $\rightarrow$ MultiQC)

```bash
pixi run bash scripts/run_all.sh \
  --indir data/raw \
  --star-index data/reference/star_index \
  --gtf data/reference/genes.gtf \
  --outdir results/bash_run1 \
  --threads 8
```

หรือใช้การตั้งค่าจาก `config.env`:
```bash
pixi run bash scripts/run_all.sh --config config.env
```

**Options ที่รองรับใน `run_all.sh`**:
| Flag | Alternative | ความหมาย | ค่า Default |
|---|---|---|---|
| `-i` | `--indir` | โฟลเดอร์เก็บไฟล์ FASTQ ดิบ | `${RAW_DIR}` (`data/raw`) |
| `-o` | `--outdir` | โฟลเดอร์ปลายทางหลักของผลลัพธ์ | Directory ปัจจุบัน หรือสร้างตามลำดับ |
| `-x` | `--star-index` | โฟลเดอร์ STAR genome index | `${STAR_INDEX}` |
| `-g` | `--gtf` | พาธไฟล์ annotation GTF | `${GTF_FILE}` |
| `-t` | `--threads` | จำนวน CPU threads | `${THREADS}` (Auto-detect) |
| `-c` | `--config` | ระบุไฟล์ config.env | `${PROJECT_ROOT}/config.env` |

> [!WARNING]
> การรัน Bash Scripts ซ้ำในโฟลเดอร์เดิมจะไม่มีระบบ Resume อัตโนมัติเหมือน Nextflow โดยจะเขียนทับไฟล์เดิม และ featureCounts อาจดึงเอาไฟล์ BAM เก่าที่ค้างอยู่มารวมด้วย ดังนั้นแนะนำให้ระบุ `--outdir` แยกต่างหากสำหรับแต่ละ run

---

### 4.3 วิธีที่ 2: รันแยกทีละขั้นตอน (Step-by-Step)

หากต้องการรันเพื่อตรวจสอบหรือแก้ไขปัญหาเฉพาะจุด สามารถรันเรียงลำดับได้ดังนี้:

#### ขั้นตอนที่ 0: FastQC ข้อมูลดิบ (Raw QC)
```bash
pixi run bash scripts/qc_data_00.sh data/raw reports/qc_raw 8
```
- **Input**: ไฟล์ FASTQ ดิบทั้งหมดใน `data/raw/`
- **Output**: รายงาน `*.html` และ `*.zip` ใน `reports/qc_raw/`

#### ขั้นตอนที่ 1: ตัด Adapter และกรองคุณภาพด้วย fastp
```bash
pixi run bash scripts/fastp_01.sh data/raw data/clean reports/fastp 8
```
- **Input**: Paired-end FASTQ (`_R1` และ `_R2`)
- **Output**: Clean reads ใน `data/clean/*_R{1,2}.clean.fastq.gz` พร้อมรายงาน JSON และ HTML ใน `reports/fastp/`

#### ขั้นตอนที่ 1.5: ตรวจสอบคุณภาพหลังตัดแต่ง (Clean QC)
```bash
pixi run bash scripts/qc_data_00.sh data/clean reports/qc_clean 8
```
- **Output**: รายงาน FastQC ของ Clean reads ใน `reports/qc_clean/` เพื่อเทียบกับ Raw QC

#### ขั้นตอนที่ 2: จัดตำแหน่งลำดับเบสด้วย STAR (Alignment)
```bash
pixi run bash scripts/align_02.sh data/clean data/aligned data/reference/star_index 8
```
- **Input**: Clean FASTQ ใน `data/clean/` + โฟลเดอร์ STAR Index
- **Output**: Coordinate-sorted BAM (`*_Aligned.sortedByCoord.out.bam`) พร้อมดัชนี `.bai` และไฟล์ Log (`*_Log.final.out`) ใน `data/aligned/`

#### ขั้นตอนที่ 3: นับจำนวน Read ต่อยีนด้วย featureCounts (Quantification)
```bash
pixi run bash scripts/quant_03.sh data/aligned data/counts data/reference/genes.gtf 8
```
- **Input**: BAM files ทั้งหมดใน `data/aligned/` + ไฟล์ `genes.gtf`
- **Output**:
  - `data/counts/gene_cell_count_matrix.tsv` (ตาราง Gene-by-Cell สำหรับ Downstream)
  - `data/counts/featurecounts_raw.txt.summary` (สถิติการ Assigned/Unassigned)

#### ขั้นตอนที่ 4: สรุปผลภาพรวมทั้งหมดด้วย MultiQC
```bash
pixi run bash scripts/multiqc_04.sh \
  reports/multiqc \
  reports/qc_raw \
  reports/qc_clean \
  reports/fastp \
  data/aligned \
  data/counts
```
- **Output**: `reports/multiqc/single_cell_multiqc_report.html`

---

## 5. Track 2: การรันด้วย Nextflow DSL2 Pipeline

Nextflow เป็นระบบ Pipeline มาตรฐานระดับ Production เหมาะสำหรับการรันข้อมูลขนาดใหญ่ มีข้อดีคือ:
- **Resume Cache**: เมื่อเกิด Error หรือต้องการรันต่อ สั่ง `-resume` จะทำงานต่อจากจุดเดิมทันทีโดยไม่ต้องเริ่มใหม่
- **Parallelization**: บริหารจัดการคิวงานและแบ่งงานตามจำนวนเซลล์อย่างมีประสิทธิภาพ
- **Container Support**: สลับไปรันบน Docker / Singularity / HPC Slurm ได้ง่าย

---

### 5.1 การเตรียม Samplesheet

Nextflow ในโปรเจกต์นี้รองรับทั้งข้อมูล **Paired-End** และ **Single-End**:

#### รูปแบบไฟล์ CSV (`samplesheet.csv`):
```csv
sample,fastq_1,fastq_2
cell_01,/data/raw/cell_01_R1.fastq.gz,/data/raw/cell_01_R2.fastq.gz
cell_02,/data/raw/cell_02_R1.fastq.gz,/data/raw/cell_02_R2.fastq.gz
cell_SE,/data/raw/cell_SE.fastq.gz,
```
- สำหรับ **Paired-End**: ใส่ทั้ง `fastq_1` และ `fastq_2`
- สำหรับ **Single-End**: ใส่เฉพาะ `fastq_1` และปล่อยคอลัมน์ `fastq_2` ให้ว่างไว้
- Path สามารถเป็น **Absolute Path** หรือ **Relative Path** (เทียบจากตำแหน่งที่วางไฟล์ CSV นั้น)
- ชื่อ `sample` ต้องไม่ซ้ำกัน และใช้เฉพาะตัวอักษรภาษาอังกฤษ ตัวเลข และเครื่องหมาย `_ . -` เท่านั้น

#### การใช้สคริปต์ช่วยสร้าง Samplesheet:
หากมีรายชื่อไฟล์ใน `docs/filenames.txt` สามารถรันสคริปต์อัตโนมัติ:
```bash
python scripts/generate_samplesheet.py docs/filenames.txt /path/to/fastq_dir samplesheet.csv
```

---

### 5.2 คำสั่งการรัน Nextflow

รันคำสั่งจาก Root directory:

```bash
pixi run nextflow run nextflow/main.nf \
  --input samplesheet.csv \
  --star_index data/reference/star_index \
  --gtf data/reference/genes.gtf \
  --outdir results/nextflow \
  --threads 8 \
  -profile local \
  -resume
```

#### กรณีต้องการรันโดยระบุ Glob Pattern แทน Samplesheet:
```bash
pixi run nextflow run nextflow/main.nf \
  --reads 'data/raw/*_R{1,2}*.fastq.gz' \
  --star_index data/reference/star_index \
  --gtf data/reference/genes.gtf \
  --outdir results/nextflow \
  --threads 8 \
  -profile local \
  -resume
```

#### ตาราง Parameters หลักของ Nextflow:
| Parameter | ความหมาย | ตัวอย่าง / ค่า Default |
|---|---|---|
| `--input` หรือ `--sample_sheet` | Path ไปยังไฟล์ samplesheet CSV | `samplesheet.csv` |
| `--reads` | File pattern กรณีไม่ใช้ samplesheet | `'data/raw/*_R{1,2}*.fastq.gz'` |
| `--star_index` | โฟลเดอร์ STAR Genome Index | `data/reference/star_index` |
| `--gtf` | Path ไปยังไฟล์ annotation GTF | `data/reference/genes.gtf` |
| `--outdir` | โฟลเดอร์เก็บผลลัพธ์ปลายทาง | `results/nextflow` |
| `--threads` | จำกัด CPU threads ต่อกระบวนการ | `4` หรือ `8` |
| `--strandedness` | ค่า Strandedness สำหรับ featureCounts | `0` (Unstranded), `1`, `2` |
| `--min_length` | ความยาวต่ำสุดหลัง trim (bp) | `35` |
| `--min_quality` | Phred quality cutoff | `20` |
| `-profile` | Execution profile | `local`, `docker`, `singularity`, `slurm` |
| `-resume` | ใช้งาน cache จากการรันครั้งก่อน | Flag |

---

### 5.3 การใช้ `-resume`

เมื่อเพิ่ม flag `-resume` ในคำสั่ง Nextflow:
1. Nextflow จะตรวจดู directory `work/` หาก input และ parameter ของ task ใดไม่เปลี่ยนแปลง Nextflow จะข้ามขั้นตอนนั้นทันที
2. หากงานหยุดชะงักจากไฟดับ หรือ Out of Memory ให้เพิ่ม RAM แล้วรันคำสั่งเดิมด้วย `-resume` ต่อได้เลย

---

### 5.4 Execution Profiles

กำหนดผ่าน `-profile <name>` ในไฟล์ `nextflow/nextflow.config`:

1. **`-profile local`** (ค่ามาตรฐานเมื่อใช้คู่กับ Pixi):
   - ใช้ binary tools ที่ติดตั้งอยู่ใน `.pixi/envs/default/bin` โดยตรง
2. **`-profile docker`**:
   - ดึง BioContainers ของแต่ละเครื่องมือมาทำงานโดยอัตโนมัติ (FastQC, fastp, STAR, Subread, MultiQC)
3. **`-profile singularity`**:
   - เหมาะกับระบบ HPC/Cluster ที่ไม่อนุญาตให้ใช้ Docker
4. **`-profile slurm`**:
   - กระจายงานเข้าคิว Slurm บน Supercomputer

---

### 5.5 การปรับแต่ง CPU/RAM ใน Nextflow

หากต้องการปรับขนาด RAM หรือ CPU สำหรับ STAR Alignment เฉพาะงาน ให้สร้างไฟล์ config เพิ่มเติม เช่น `custom_resources.config`:

```groovy
process {
    withName: 'STAR_ALIGN' {
        cpus   = 8
        memory = 48.GB
        maxForks = 2 // รัน STAR พร้อมกันได้สูงสุด 2 เซลล์เพื่อป้องกัน RAM หมด
    }
}
params.star_sort_ram = 8000000000 // 8 GB สำหรับ samtools coordinate sort
```

จากนั้นรัน Nextflow โดยเพิ่ม `-c custom_resources.config`:
```bash
pixi run nextflow run nextflow/main.nf \
  -c custom_resources.config \
  --input samplesheet.csv \
  --star_index data/reference/star_index \
  --gtf data/reference/genes.gtf \
  --outdir results/nextflow \
  -profile local -resume
```

---

## 6. โครงสร้างผลลัพธ์และการตรวจสอบคุณภาพ

### 6.1 โครงสร้างโฟลเดอร์ผลลัพธ์ (Output Tree)

```text
results/nextflow/ (หรือ results/bash/)
├── counts/
│   ├── gene_cell_count_matrix.tsv           # ⭐ ผลลัพธ์หลัก: เมทริกซ์ยีน x เซลล์ (Tab-separated)
│   ├── featurecounts_raw.txt                # ข้อมูลดิบจาก featureCounts
│   └── featurecounts_raw.txt.summary        # สรุปจำนวน read ที่ถูกนับ / ไม่ถูกนับ
├── aligned/
│   ├── <sample>_Aligned.sortedByCoord.out.bam # BAM ไฟล์ที่เรียงลำดับพิกัดแล้ว
│   ├── <sample>_Aligned.sortedByCoord.out.bam.bai # Index ของ BAM
│   └── <sample>_Log.final.out               # สถิติ Alignment จาก STAR
├── fastp/
│   ├── <sample>_fastp.html                  # รายงานกราฟคุณภาพและ adapter trimming
│   └── <sample>_fastp.json                  # สถิติแบบ JSON
├── fastqc/
│   ├── raw/                                 # FastQC ก่อน trim
│   └── clean/                               # FastQC หลัง trim
├── multiqc/
│   └── single_cell_multiqc_report.html      # ⭐ รายงานสรุปคุณภาพของทุกเซลล์ในหน้าเดียว
└── pipeline_info/                           # (เฉพาะ Nextflow)
    ├── execution_report.html                # รายงานการใช้ CPU/RAM ของแต่ละ Task
    ├── execution_timeline.html              # ไทม์ไลน์เวลาที่ใช้
    └── execution_trace.txt                  # รายละเอียด Exit Code และ Duration
```

### 6.2 จุดที่ต้องตรวจสอบใน MultiQC Report

เปิดไฟล์ `single_cell_multiqc_report.html` ด้วยเว็บเบราว์เซอร์:
1. **fastp / Adapter Trimming**:
   - อัตรา % Passed Filter ควรสูงกว่า 80–90%
   - ไม่ควรมี Adapter contamination หลงเหลือใน Clean FastQC
2. **STAR Alignment**:
   - **Uniquely mapped reads %**: สำหรับ Smart-seq2 คุณภาพดี ควรอยู่ที่ **> 60–80%**
   - หาก **% Unmapped too short** สูงผิดปกติ อาจเกิดจาก Degradation ของ RNA หรือ Read สั้นเกินไป
   - หาก **% Unmapped other** สูง อาจเกิดจากการเลือก Reference Genome ผิดสปีชีส์ หรือมี DNA การปนเปื้อน
3. **featureCounts Assignment**:
   - อัตรา **Assigned reads** ต่อยีนควรอยู่ระหว่าง **40–70%** (ส่วนที่เหลือมักเป็น Unassigned_NoFeatures คือ Intron/Intergenic หรือ Unassigned_Ambiguity)

---

## 7. การส่งต่อไปยัง Downstream Analysis

หลังจากได้ไฟล์ `gene_cell_count_matrix.tsv` สามารถเปิดใช้งาน JupyterLab ในโมดูล `workbench` เพื่อวิเคราะห์ข้อมูล Single-Cell ด้วย AnnData / Scanpy / PCA / UMAP:

```bash
cd workbench
pixi run lab
```

เปิดสมุดบันทึก `workbench/notebooks/00_end_to_end.ipynb` แล้วเปลี่ยนพาธตัวแปร:
```python
# ชี้ไปยังผลลัพธ์ของ Upstream Pipeline
COUNTS_PATH = "../results/nextflow/counts/gene_cell_count_matrix.tsv"
# หรือสำหรับ Bash:
# COUNTS_PATH = "../results/bash/data/counts/gene_cell_count_matrix.tsv"
```

---

## 8. ปัญหาที่พบบ่อยและวิธีแก้ไข (Troubleshooting & FAQs)

### 🔴 1. STAR ขึ้น Error "Exit Code 137" หรือ "Killed"
- **สาเหตุ**: หน่วยความจำของเครื่อง (RAM) ไม่เพียงพอ การโหลด Human STAR Index เข้าหน่วยความจำต้องใช้ RAM อย่างน้อย ~32 GB
- **วิธีแก้ไข**:
  1. ตรวจสอบ RAM คงเหลือด้วย `free -h`
  2. หากรัน Nextflow ให้จำกัด `maxForks = 1` สำหรับ process `STAR_ALIGN` เพื่อไม่ให้รัน STAR พร้อมกันหลายเซลล์
  3. ปรับลด `params.star_sort_ram` ให้เหมาะสม (เช่น `4000000000` = 4GB)

### 🔴 2. FastQC ขึ้น "OutOfMemoryError: Java heap space"
- **สาเหตุ**: Java heap memory สำหรับ FastQC มีไม่พอ
- **วิธีแก้ไข**: รันคำสั่งเพิ่ม Heap ก่อนเริ่มรันสคริปต์:
  ```bash
  export _JAVA_OPTIONS="-Xmx4096m"
  ```

### 🔴 3. fastp ขึ้น "R1 and R2 have different number of reads"
- **สาเหตุ**: ไฟล์ Fastq คู่ R1 และ R2 ไม่สมบูรณ์ หรือดาวน์โหลดมาไม่ครบ
- **วิธีแก้ไข**: ตรวจสอบความสมบูรณ์ของ gzip:
  ```bash
  gzip -t sample_R1.fastq.gz
  gzip -t sample_R2.fastq.gz
  ```

### 🔴 4. featureCounts ได้ Assigned reads เป็น 0% หรือต่ำมาก (< 10%)
- **สาเหตุ**:
  1. ค่า **Strandedness** ไม่ถูกต้อง (Smart-seq2 ส่วนใหญ่เป็น unstranded = `0`) หากเผลอตั้งเป็น `1` หรือ `2` โปรแกรมจะนับผิด strand
  2. รูปแบบ Chromosome naming ใน GTF ไม่ตรงกับ STAR Index (เช่น ฝั่งหนึ่งเป็น `chr1`, อีกฝั่งเป็น `1`)
- **วิธีแก้ไข**: ตรวจสอบ header ของ BAM และ GTF:
  ```bash
  samtools view -H sample_Aligned.sortedByCoord.out.bam | grep '@SQ' | head -n 3
  grep -v '^#' data/reference/genes.gtf | cut -f1 | head -n 3
  ```
  หากพบว่าข้างหนึ่งมี `chr` แต่อีกข้างไม่มี ให้ใช้ GTF และ FASTA จากแหล่งเดียวกัน (เช่น Ensembl เหมือนกันทั้งคู่ หรือ GENCODE เหมือนกันทั้งคู่)

### 🔴 5. Nextflow ฟ้อง "Incomplete STAR index: missing Genome/SA/SAindex"
- **สาเหตุ**: ในโฟลเดอร์ STAR Index ขาดไฟล์สำคัญอันใดอันหนึ่งเนื่องจากสร้าง index ไม่สำเร็จ
- **วิธีแก้ไข**: ตรวจสอบโฟลเดอร์ index ต้องมีไฟล์ครบทั้ง 3 ไฟล์:
  ```bash
  ls -lh data/reference/star_index/Genome data/reference/star_index/SA data/reference/star_index/SAindex
  ```
  หากไม่มี ให้สั่งรัน `scripts/setup_reference.sh` ใหม่อีกครั้งให้เสร็จสมบูรณ์
