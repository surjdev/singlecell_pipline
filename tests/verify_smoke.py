"""Check both real-tool smoke runs: python tests/verify_smoke.py /tmp/sc-demo."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('outdir', type=Path)
root = parser.parse_args().outdir
expected = 'Geneid\tcell_a\tcell_b\ngene1\t200\t200\n'
for run, counts, report in [
    ('bash-results', 'data/counts', 'reports/multiqc'),
    ('nf-results', 'counts', 'multiqc'),
]:
    base = root / run
    assert (base / counts / 'gene_cell_count_matrix.tsv').read_text() == expected
    assert (base / report / 'single_cell_multiqc_report.html').stat().st_size > 1000
    qc = base / report / 'single_cell_multiqc_report_data/multiqc_fastqc.txt'
    assert len(qc.read_text().splitlines()) == 9, 'Expected 8 raw/clean mate reports'
print('Both routes produced expected counts and complete raw/clean QC reports')
