"""
Generates the 3 interactive Jupyter Notebooks for the single-cell pipeline.
"""

import nbformat as nbf
from pathlib import Path

def create_notebook_01():
    nb = nbf.v4.new_notebook()
    cells = []
    
    # Header
    cells.append(nbf.v4.new_markdown_cell("""# 🔬 01. Single-Cell RAW FASTQ Quality Control & Exploration
### Single-Cell RNA-Seq (10x Genomics 3' Chromium)
This notebook performs quality control inspection of RAW single-cell sequencing reads prior to preprocessing.

#### Key Objectives:
1. Understand 10x Chromium read architecture (R1: 16bp Cell Barcode + 12bp UMI; R2: cDNA).
2. Evaluate per-cycle Phred quality scores across Cell Barcode, UMI, and cDNA.
3. Compute Cell Barcode diversity and render the **Knee Plot** (Rank vs Read Depth).
4. Measure Whitelist match rates and identify error-correctable barcode populations.
"""))

    # Imports & Setup
    cells.append(nbf.v4.new_code_cell("""import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path("..").resolve()))

import json
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set plot styles
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['font.size'] = 11

from src.utils import load_config
from src.sc_qc import SingleCellQC

# Load configuration
config = load_config("../config/pipeline_config.yaml")
paths = config["paths"]
print(f"Loaded config for project: {config['project']['name']}")
"""))

    # SingleCellQC run
    cells.append(nbf.v4.new_markdown_cell("""## 1. Run Single-Cell QC Analysis
Let's initialize the `SingleCellQC` engine on the raw FASTQ files.
"""))

    cells.append(nbf.v4.new_code_cell("""qc_raw = SingleCellQC(
    r1_path=f"../{paths['sample_r1']}",
    r2_path=f"../{paths['sample_r2']}",
    cb_len=config['chemistry']['r1_structure']['cell_barcode_len'],
    umi_len=config['chemistry']['r1_structure']['umi_len'],
    whitelist_path=f"../{paths['whitelist_file']}",
    sample_id="raw_sample_01"
)

raw_results = qc_raw.analyze_fastq()
metrics = raw_results["summary_metrics"]
pd.DataFrame([metrics]).T.rename(columns={0: "Raw Value"})
"""))

    # Plot 1: Per-cycle quality
    cells.append(nbf.v4.new_markdown_cell("""## 2. Per-Cycle Sequencing Quality Profile
In single-cell 10x sequencing:
- **Cycles 1-16 (R1):** Cell Barcode
- **Cycles 17-28 (R1):** Unique Molecular Identifier (UMI)
- **Cycles 1-91 (R2):** cDNA Transcript read
"""))

    cells.append(nbf.v4.new_code_cell("""fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.5))

# R1 Profile
r1_quals = raw_results["r1_mean_qual_per_cycle"]
ax1.plot(range(1, len(r1_quals) + 1), r1_quals, color="#1f77b4", lw=2, marker='o', markersize=4)
ax1.axvspan(1, 16, color="#4CAF50", alpha=0.15, label="Cell Barcode (1-16bp)")
ax1.axvspan(16, 28, color="#FF9800", alpha=0.15, label="UMI (17-28bp)")
ax1.axhline(30, color="red", linestyle="--", alpha=0.7, label="Q30 Threshold")
ax1.set_title("R1: Cell Barcode & UMI Quality", fontweight='bold')
ax1.set_xlabel("Sequencing Cycle (bp)")
ax1.set_ylabel("Mean Phred Quality Score")
ax1.set_ylim(0, 42)
ax1.legend(loc="lower right")

# R2 Profile
r2_quals = raw_results["r2_mean_qual_per_cycle"]
ax2.plot(range(1, len(r2_quals) + 1), r2_quals, color="#e91e63", lw=2)
ax2.axhline(30, color="red", linestyle="--", alpha=0.7, label="Q30 Threshold")
ax2.set_title("R2: cDNA Quality (3' Degradation)", fontweight='bold')
ax2.set_xlabel("Sequencing Cycle (bp)")
ax2.set_ylabel("Mean Phred Quality Score")
ax2.set_ylim(0, 42)
ax2.legend(loc="lower left")

plt.tight_layout()
plt.show()
"""))

    # Plot 2: Knee Plot
    cells.append(nbf.v4.new_markdown_cell("""## 3. Cell Barcode Knee Plot (Barcode Rank vs Read Depth)
The Knee Plot is the foundational QC plot in single-cell genomics:
- High-depth plateau represents **authentic single cells**.
- The steep drop (inflection / knee point) separates captured single cells from **ambient background RNA droplets**.
"""))

    cells.append(nbf.v4.new_code_cell("""counts = np.array(raw_results["barcode_rank_counts"])
ranks = np.arange(1, len(counts) + 1)
knee_stats = raw_results["knee_stats"]
est_cells = knee_stats["estimated_cells"]

plt.figure(figsize=(8, 5.5))
plt.loglog(ranks, counts, color="#3f51b5", lw=2.5, label="Barcode Depth Curve")
plt.axvline(est_cells, color="#e53935", linestyle="--", lw=2, 
            label=f"Knee Cutoff: ~{est_cells} cells ({knee_stats['fraction_reads_in_cells']*100:.1f}% reads)")
plt.scatter([est_cells], [counts[min(est_cells-1, len(counts)-1)]], color="red", s=80, zorder=5)

plt.title("Single-Cell Barcode Rank vs. UMI/Read Count (Knee Plot)", fontweight='bold')
plt.xlabel("Barcode Rank (log10)")
plt.ylabel("Read Depth per Barcode (log10)")
plt.legend(frameon=True, facecolor="white", loc="lower left")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()
"""))

    # Whitelist match rate
    cells.append(nbf.v4.new_markdown_cell("""## 4. Barcode Whitelist Match Rate Breakdown
Let's see what percentage of raw cell barcodes match the official whitelist exactly, require 1-bp error correction, or are invalid chimeric reads.
"""))

    cells.append(nbf.v4.new_code_cell("""wl_labels = ["Exact Whitelist Match", "1-bp Mismatch (Recoverable)", "Invalid / Noise"]
wl_values = [
    metrics.get("whitelist_exact_match_rate", 0) * 100,
    metrics.get("whitelist_1bp_mismatch_rate", 0) * 100,
    metrics.get("whitelist_invalid_rate", 0) * 100
]
colors = ["#4CAF50", "#FF9800", "#F44336"]

fig, ax = plt.subplots(figsize=(7, 5))
wedges, texts, autotexts = ax.pie(
    wl_values, labels=wl_labels, autopct='%1.2f%%',
    colors=colors, startangle=140, explode=(0, 0.05, 0.1),
    wedgeprops=dict(width=0.6, edgecolor='w')
)
plt.setp(autotexts, size=10, weight="bold")
plt.title("Raw Cell Barcode Whitelist Conformity", fontweight='bold')
plt.tight_layout()
plt.show()
"""))

    nb.cells = cells
    return nb

def create_notebook_02():
    nb = nbf.v4.new_notebook()
    cells = []
    
    # Header
    cells.append(nbf.v4.new_markdown_cell("""# ⚙️ 02. Single-Cell Preprocessing & Barcode Error Correction
### Trimming, Poly-G / Poly-A Filtering, and 1-Hamming Distance Error Correction

#### Pipeline Steps:
1. **Quality & Adapter Trimming with `fastp`**: Sliding window Phred filtering, removal of NextSeq/NovaSeq Poly-G dark cycle artifacts, and TSO adapter removal.
2. **Whitelist Error Correction**: 1-Hamming distance error-correction of 1-bp sequencing mutations in cell barcodes.
3. **Clean FASTQ Generation**: Outputting standard paired FASTQ and extracted FASTQ format (`@READ_CB_UMI`).
"""))

    # Setup
    cells.append(nbf.v4.new_code_cell("""import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path("..").resolve()))

import json
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from src.utils import load_config
from src.sc_preprocess import SingleCellPreprocessor

config = load_config("../config/pipeline_config.yaml")
print(f"Loaded configuration for sample: {config['project']['sample_id']}")
"""))

    # Fastp Trimming Execution
    cells.append(nbf.v4.new_markdown_cell("""## 1. Execute `fastp` Quality & Adapter Trimming
We run `fastp` with single-cell parameters:
- Preserve R1 (28bp CB+UMI)
- Trim low-quality bases from R2 3' end
- Strip poly-G and poly-A tails
"""))

    cells.append(nbf.v4.new_code_cell("""preprocessor = SingleCellPreprocessor(config)

temp_r1 = Path("../data/temp/nb_trimmed_R1.fastq.gz")
temp_r2 = Path("../data/temp/nb_trimmed_R2.fastq.gz")
fp_html = Path("../reports/fastp/nb_fastp.html")
fp_json = Path("../reports/fastp/nb_fastp.json")

fp_stats = preprocessor.run_fastp_trimming(
    r1_in=f"../{config['paths']['sample_r1']}",
    r2_in=f"../{config['paths']['sample_r2']}",
    r1_out=temp_r1,
    r2_out=temp_r2,
    report_html=fp_html,
    report_json=fp_json
)

print(f"Fastp filtering completed!")
print(f"Reads before: {fp_stats['summary']['before_filtering']['total_reads']:,}")
print(f"Reads after:  {fp_stats['summary']['after_filtering']['total_reads']:,}")
print(f"Q30 rate before: {fp_stats['summary']['before_filtering']['q30_rate']*100:.2f}%")
print(f"Q30 rate after:  {fp_stats['summary']['after_filtering']['q30_rate']*100:.2f}%")
"""))

    # Poly-G / Poly-A Trimming Breakdown
    cells.append(nbf.v4.new_markdown_cell("""## 2. Poly-X & Adapter Trimming Results
Let's inspect the exact number of reads trimmed for Poly-A / Poly-G artifacts.
"""))

    cells.append(nbf.v4.new_code_cell("""polyx_info = fp_stats.get("polyx_trimming", {})
polyx_reads = polyx_info.get("polyx_trimmed_reads", {})

df_polyx = pd.DataFrame(list(polyx_reads.items()), columns=["Base", "Trimmed Reads"]).sort_values("Trimmed Reads", ascending=False)

plt.figure(figsize=(7, 4))
sns.barplot(data=df_polyx, x="Base", y="Trimmed Reads", palette="viridis")
plt.title("Reads Trimmed by Base Type (Poly-X Filter)", fontweight='bold')
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.show()
"""))

    # Barcode Error Correction
    cells.append(nbf.v4.new_markdown_cell("""## 3. Barcode Error Correction (1-Hamming Distance)
We now map each cell barcode against the 10x whitelist:
- If exact match $\\rightarrow$ Keep as-is.
- If 1-bp mismatch with unambiguous whitelist barcode $\\rightarrow$ Error-correct to true barcode.
- If $>1$-bp mismatch or ambiguous $\\rightarrow$ Discard.
"""))

    cells.append(nbf.v4.new_code_cell("""clean_r1 = Path(f"../{config['paths']['clean_r1']}")
clean_r2 = Path(f"../{config['paths']['clean_r2']}")
extracted_out = Path(f"../{config['paths']['extracted_fastq']}")

bc_stats = preprocessor.filter_and_correct_barcodes(
    r1_in=temp_r1,
    r2_in=temp_r2,
    r1_out=clean_r1,
    r2_out=clean_r2,
    extracted_out=extracted_out,
    whitelist_file=f"../{config['paths']['whitelist_file']}"
)

# Clean up temporary intermediate file
if temp_r1.exists(): temp_r1.unlink()
if temp_r2.exists(): temp_r2.unlink()

pd.DataFrame([bc_stats]).T.rename(columns={0: "Count / Rate"})
"""))

    # Barcode Yield Visualization
    cells.append(nbf.v4.new_markdown_cell("""## 4. Preprocessing Read Yield Breakdown"""))

    cells.append(nbf.v4.new_code_cell("""categories = ["Exact Whitelist", "1-bp Error-Corrected", "Discarded Invalid"]
vals = [
    bc_stats["exact_whitelist_matches"],
    bc_stats["corrected_1bp_mismatches"],
    bc_stats["discarded_invalid_barcodes"]
]

plt.figure(figsize=(8, 4.5))
bars = plt.bar(categories, vals, color=["#2e7d32", "#1976d2", "#d32f2f"], width=0.55)
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + max(vals)*0.02, f"{yval:,} ({yval/bc_stats['input_reads']*100:.1f}%)", ha='center', fontweight='bold')

plt.title("Barcode Recovery & Filtering Summary", fontweight='bold')
plt.ylabel("Read Count")
plt.ylim(0, max(vals) * 1.15)
plt.tight_layout()
plt.show()
"""))

    nb.cells = cells
    return nb

def create_notebook_03():
    nb = nbf.v4.new_notebook()
    cells = []
    
    # Header
    cells.append(nbf.v4.new_markdown_cell("""# 📊 03. Post-QC, MultiQC Dashboard & Benchmark Analysis
### Comparative Analysis: RAW FASTQ vs Clean FASTQ

#### Key Objectives:
1. Quantify quality gains ($Q30$, Poly-A/G removal, adapter cleanup).
2. Compare Knee Plot curves before and after barcode error-correction.
3. Review interactive **MultiQC** dashboard.
4. Verify downstream alignment readiness for **STARsolo / Kallisto-Bustools / CellRanger**.
"""))

    # Imports & Setup
    cells.append(nbf.v4.new_code_cell("""import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path("..").resolve()))

import json
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import IFrame, display, HTML

from src.utils import load_config
from src.sc_qc import SingleCellQC

config = load_config("../config/pipeline_config.yaml")
paths = config["paths"]
"""))

    # Run Clean QC
    cells.append(nbf.v4.new_markdown_cell("""## 1. Compute Post-Clean QC Metrics"""))

    cells.append(nbf.v4.new_code_cell("""# Run QC on Clean FASTQ
qc_clean = SingleCellQC(
    r1_path=f"../{paths['clean_r1']}",
    r2_path=f"../{paths['clean_r2']}",
    cb_len=config['chemistry']['r1_structure']['cell_barcode_len'],
    umi_len=config['chemistry']['r1_structure']['umi_len'],
    whitelist_path=f"../{paths['whitelist_file']}",
    sample_id="clean_sample_01"
)
clean_results = qc_clean.analyze_fastq()

# Load raw QC report
with open(f"../{paths['qc_raw_dir']}/{config['project']['sample_id']}_raw_sc_qc.json", "r") as f:
    raw_saved = json.load(f)
raw_m = raw_saved["summary_metrics"]
clean_m = clean_results["summary_metrics"]
"""))

    # Comparison Table
    cells.append(nbf.v4.new_markdown_cell("""## 2. Before vs After Quality Benchmark Table"""))

    cells.append(nbf.v4.new_code_cell("""comparison_df = pd.DataFrame([
    {"Metric": "Total Reads", "RAW FASTQ": f"{raw_m['total_reads']:,}", "CLEAN FASTQ": f"{clean_m['total_reads']:,}", "Change": f"{(clean_m['total_reads']-raw_m['total_reads'])/raw_m['total_reads']*100:.1f}%"},
    {"Metric": "cDNA (R2) Q30 Rate", "RAW FASTQ": f"{raw_m['r2_q30_fraction']*100:.2f}%", "CLEAN FASTQ": f"{clean_m['r2_q30_fraction']*100:.2f}%", "Change": f"+{(clean_m['r2_q30_fraction']-raw_m['r2_q30_fraction'])*100:.2f}%"},
    {"Metric": "Barcode Q30 Rate", "RAW FASTQ": f"{raw_m['cb_q30_fraction']*100:.2f}%", "CLEAN FASTQ": f"{clean_m['cb_q30_fraction']*100:.2f}%", "Change": f"{(clean_m['cb_q30_fraction']-raw_m['cb_q30_fraction'])*100:.2f}%"},
    {"Metric": "Valid Barcodes", "RAW FASTQ": f"{raw_m.get('valid_barcode_fraction',0)*100:.2f}%", "CLEAN FASTQ": f"{clean_m.get('valid_barcode_fraction',0)*100:.2f}%", "Change": "100% Whitelist"},
    {"Metric": "Poly-A Artifact Rate", "RAW FASTQ": f"{raw_m['poly_a_rate']*100:.2f}%", "CLEAN FASTQ": f"{clean_m['poly_a_rate']*100:.2f}%", "Change": "Eliminated"},
    {"Metric": "Poly-G Artifact Rate", "RAW FASTQ": f"{raw_m['poly_g_rate']*100:.2f}%", "CLEAN FASTQ": f"{clean_m['poly_g_rate']*100:.2f}%", "Change": "Eliminated"}
])
display(comparison_df.style.set_properties(**{'text-align': 'left'}))
"""))

    # Comparison Plots
    cells.append(nbf.v4.new_markdown_cell("""## 3. Comparative Visualizations"""))

    cells.append(nbf.v4.new_code_cell("""fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Q30 Comparison Barplot
metrics_names = ["cDNA Q30", "CB Q30", "UMI Q30"]
raw_vals = [raw_m['r2_q30_fraction']*100, raw_m['cb_q30_fraction']*100, raw_m['umi_q30_fraction']*100]
clean_vals = [clean_m['r2_q30_fraction']*100, clean_m['cb_q30_fraction']*100, clean_m['umi_q30_fraction']*100]

x = np.arange(len(metrics_names))
width = 0.35

ax1.bar(x - width/2, raw_vals, width, label='Raw FASTQ', color='#90a4ae')
ax1.bar(x + width/2, clean_vals, width, label='Clean FASTQ', color='#2e7d32')
ax1.set_ylabel('Q30 Percentage (%)')
ax1.set_title('Phred Quality (Q30) Before vs After Cleaning', fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(metrics_names)
ax1.set_ylim(0, 110)
ax1.legend()

# Knee Plot Comparison
raw_counts = np.array(raw_saved["top_100_barcode_counts"])
clean_counts = np.array(clean_results["barcode_rank_counts"][:100])

ax2.plot(range(1, len(raw_counts)+1), raw_counts, label='Raw (with uncorrected)', color='#f57c00', lw=2)
ax2.plot(range(1, len(clean_counts)+1), clean_counts, label='Clean (Recovered & Corrected)', color='#1976d2', lw=2)
ax2.set_title('Knee Plot Comparison (Top Barcodes)', fontweight='bold')
ax2.set_xlabel('Barcode Rank')
ax2.set_ylabel('Read Depth')
ax2.legend()

plt.tight_layout()
plt.show()
"""))

    # MultiQC Dashboard Link
    cells.append(nbf.v4.new_markdown_cell("""## 4. MultiQC Aggregated Quality Dashboard
MultiQC compiles FastQC, Fastp, and tool metrics into an interactive HTML report.
"""))

    cells.append(nbf.v4.new_code_cell("""multiqc_path = Path("../reports/multiqc/single_cell_multiqc_report.html")
if multiqc_path.exists():
    print(f"✔ MultiQC report ready at: {multiqc_path.resolve()}")
    display(HTML(f'<p>👉 Open the interactive report: <a href="../{paths[\"multiqc_dir\"]}/single_cell_multiqc_report.html" target="_blank" style="color: #1976d2; font-weight: bold;">MultiQC HTML Report</a></p>'))
else:
    print("MultiQC report not found. Run `pixi run multiqc-report` to generate it.")
"""))

    # Downstream Alignment Ready Check
    cells.append(nbf.v4.new_markdown_cell("""## 5. Downstream Alignment Compatibility
The output FASTQ files are ready for standard single-cell aligners:
- **STARsolo**:
  ```bash
  STAR --genomeDir /path/to/ref --readFilesIn data/clean/scRNA_sample_01_val_R2.fastq.gz data/clean/scRNA_sample_01_val_R1.fastq.gz --soloType CB_UMI_Simple --soloCBwhitelist data/whitelist/737K-august-2016.txt
  ```
- **Kallisto / Bustools**:
  ```bash
  kallisto bus -i transcripts.idx -o bus_output/ -x 10xv3 -t 4 data/clean/scRNA_sample_01_val_R1.fastq.gz data/clean/scRNA_sample_01_val_R2.fastq.gz
  ```
- **Alevin / Salmon / UMI-tools**:
  Compatible with extracted header format in `data/clean/scRNA_sample_01_extracted.fastq.gz`.
"""))

    nb.cells = cells
    return nb

def main():
    nb_dir = Path("notebooks")
    nb_dir.mkdir(parents=True, exist_ok=True)
    
    nb1 = create_notebook_01()
    with open(nb_dir / "01_raw_data_qc.ipynb", "w", encoding="utf-8") as f:
        nbf.write(nb1, f)
    print("✔ Created notebooks/01_raw_data_qc.ipynb")
    
    nb2 = create_notebook_02()
    with open(nb_dir / "02_preprocessing_pipeline.ipynb", "w", encoding="utf-8") as f:
        nbf.write(nb2, f)
    print("✔ Created notebooks/02_preprocessing_pipeline.ipynb")
    
    nb3 = create_notebook_03()
    with open(nb_dir / "03_post_qc_and_benchmarking.ipynb", "w", encoding="utf-8") as f:
        nbf.write(nb3, f)
    print("✔ Created notebooks/03_post_qc_and_benchmarking.ipynb")

if __name__ == "__main__":
    main()
