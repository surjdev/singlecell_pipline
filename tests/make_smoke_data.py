from pathlib import Path
import random,gzip
import argparse
parser = argparse.ArgumentParser(description='Generate synthetic paired reads and a tiny reference for software smoke tests')
parser.add_argument('outdir', type=Path)
p = parser.parse_args().outdir
p.mkdir(parents=True, exist_ok=True)
r=random.Random(12); seq=''.join(r.choices('ACGT',k=12000))
(p/'genome.fa').write_text('>chr1\n'+seq+'\n')
(p/'genes.gtf').write_text('chr1\tsynthetic\texon\t1001\t11000\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "DEMO1";\n')
for cell in ['cell_a','cell_b']:
 with gzip.open(p/f'{cell}_R1.fastq.gz','wt') as a,gzip.open(p/f'{cell}_R2.fastq.gz','wt') as b:
  for i in range(200):
   j=r.randint(1100,10500); x=seq[j:j+75]; y=seq[j+125:j+200].translate(str.maketrans('ACGT','TGCA'))[::-1]
   a.write(f'@read{i}/1\n{x}\n+\n'+ 'I'*75+'\n'); b.write(f'@read{i}/2\n{y}\n+\n'+'I'*75+'\n')
(p/'samples.csv').write_text('sample,fastq_1,fastq_2\n'+''.join(f'{c},{c}_R1.fastq.gz,{c}_R2.fastq.gz\n' for c in ['cell_a','cell_b']))
(p/'index').mkdir(exist_ok=True)
