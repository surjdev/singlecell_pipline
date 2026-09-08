# Code quality and usability review — 2026-09-07

ตรวจจาก source, regression tests และการรันเครื่องมือจริงบน synthetic data
แก้ต่อจากไฟล์ที่มีการปรับอยู่ก่อนแล้ว โดยไม่ย้อนงานเดิมใน PLAN.md

| ปัญหาที่พบ | การปรับ |
|---|---|
| Bash `--outdir` ไม่เปลี่ยน derived paths | resolve paths ก่อนรันและ export config ให้ child stages |
| config/CLI override ไม่แน่นอน | config โหลดก่อน defaults; CLI ชนะ และ child ไม่โหลดซ้ำ |
| reference script มี EOF เกินและ download อาจค้างเป็นไฟล์ดูเหมือนสำเร็จ | ลบคำสั่งเกิน, curl fail/retry, decompress ลง partial ก่อน rename |
| FastQC xargs แยก path ที่มีช่องว่าง | ใช้ NUL separator |
| symlink FASTQ ไม่ถูกค้นพบ | ใช้ find -L |
| fastp sample ซ้ำอาจเขียนทับ | ตรวจ duplicate ID; cap threads ที่ 16 |
| ไม่มี clean FastQC ใน Bash runner | เพิ่มหลัง fastp และรวมใน MultiQC |
| featureCounts header ผูกกับเลขบรรทัด 2 | ตรวจ Geneid แทนเลขบรรทัด |
| Nextflow ไม่มี samplesheet validation | ตรวจ IDs/fields/mates/duplicate paths ก่อน dispatch |
| ลำดับ BAM ขึ้นกับ task completion | sort ก่อนรวม counts |
| Nextflow FastQC ชื่อ input ซ้ำชนกันตอน publish | เติม sample ID และ raw/clean ใน report names |
| MultiQC staging อาจชน basename | แยก input subdirectories และรักษา full sample names |
| Nextflow CPU option ไม่ถูกใช้ | ใช้ threads cap และให้ตัวอย่าง resource override |
| loader transpose=False ยัง transpose | ทำ orientation ให้ตรง API และ validate matrix/IDs |
| metadata ผูกด้วยตำแหน่งแทน ID | align Pandas index พร้อมตรวจ missing/duplicate |
| QC หลัง normalize ใช้ค่าที่เปลี่ยนแล้ว | ใช้ counts layer; gene_name รองรับ symbol annotation |
| normalize ซ้ำใช้ log expression | เริ่มจาก counts layer ทุกครั้ง |
| PCA ไม่จำกัดตามจำนวน HVG | จำกัด components ตาม cells และ selected genes |
| seurat_v3 HVG รับ log expression | ส่ง counts layer ตาม Scanpy contract |
| layer/representation พิมพ์ผิดกลับใช้ข้อมูลอื่น | raise error ชัดเจน |
| sklearn แปลง sparse เป็น dense โดยไม่จำเป็น | ส่ง matrix ให้ estimator โดยตรง |
| classifier report เป็น training accuracy | ใช้ out-of-fold predictions; documented preprocessing leakage caveat |
| Scrublet/Enrichr failure ถูกแสดงเป็นผลว่างหรือไม่มี doublet | raise exception แทนสร้างผลสำเร็จเทียม |
| pseudotime เลือก root เองโดยไม่มีบริบท | ต้องระบุ root cell/cluster |
| README เรียก functions ที่ไม่มีจริงและอ้าง production เกินหลักฐาน | เปลี่ยน quickstart เป็น API จริง; แยกข้อจำกัด/ผลทดสอบ |

## Verification

- Python tests: 17 tests passed; existing suite plus orientation, metadata alignment, invalid values, counts preservation,
  repeated normalization, missing layers/representations and small HVG PCA regression tests.
- Real Bash pipeline: synthetic paired FASTQ 2 cells ผ่าน FastQC/fastp/STAR/samtools/featureCounts/MultiQC.
- Real Nextflow local workflow: same input/reference, count matrix ตรงกับ Bash (gene1 = 200 ต่อเซลล์).
- New notebook: execute ทุก cell ด้วย Jupyter kernel; synthetic 120 cells × 300 genes,
  รวม markers, Pandas/NumPy/Scikit-learn และ H5AD reload assertion.
- See tests/README.md for reproduction commands.

## Remaining limits

ยังไม่ได้ทดสอบ full human/mouse FASTQ, reference download หลาย GB, Slurm, Docker หรือ Singularity
container tags และ scheduler config ต้องตรวจใน environment ปลายทางก่อนใช้
Bash ไม่รองรับ resume หรือ lane merging; ใช้ output ใหม่ต่อ dataset
TSV ingestion และ CSV/dense exports ยังใช้ RAM ตาม matrix ขนาดเต็ม: สำหรับข้อมูลใหญ่ควรใช้ MTX/H5AD/sparse
notebooks 01–05 เป็นตัวอย่างเดิมที่ต้องตั้ง input/parameters; 00 เป็น entry point ที่ทดสอบ end-to-end
Harmony, external enrichment และ biological trajectory ไม่ได้มี validation บนข้อมูลจริงในรอบนี้
QC thresholds และ cell types ต้องใช้ protocol/tissue/replicate metadata ไม่อนุมานจาก cluster ID

References: [Scanpy HVG](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.highly_variable_genes.html),
[Nextflow process inputs](https://nextflow.io/docs/stable/process.html).
