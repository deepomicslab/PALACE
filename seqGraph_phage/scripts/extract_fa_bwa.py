import subprocess
import os
import re
import sys
# Define the files and the working directory
file_a = sys.argv[1]
file_b = sys.argv[2]
fq1 = sys.argv[3]
fq2 = sys.argv[4]
work_dir = sys.argv[5]
samtools = sys.argv[6]
aligner = sys.argv[7]
tgs = int(sys.argv[8])

# Read file A and create the string of regions to extract
line_all = set()
with open(file_a, 'r') as f:
    for line in f:
        line_arr = re.split('[-+\t]', line.strip())
        line_all.update(line_arr[:-1])
        #line_all.add(' '.join(line_arr[:-1]))
line_str = ' '.join(line_all)

# Extract the subfasta
cmd = f'{samtools} faidx {file_b} {line_str} > {work_dir}/need_second.fasta'
subprocess.run(cmd, shell=True, check=True)

# Index the subfasta
if tgs != 1:
    cmd = f'{aligner} index {work_dir}/need_second.fasta'
    subprocess.run(cmd, shell=True, check=True)

# Align the fastq files to the subfasta
if not os.path.exists(f'{work_dir}/need_second.sam'):
    if tgs == 0:
        cmd = f'{aligner} mem -t 8 {work_dir}/need_second.fasta {fq1} {fq2} > {work_dir}/need_second.sam'
    else:
        cmd = f'{aligner} -x map-ont -t 64 -a -p 0.5 -N 10 --sam-hit-only -L -z 1000 -Q {work_dir}/need_second.fasta {fq1} > {work_dir}/need_second.sam'
    subprocess.run(cmd, shell=True, check=True)

# Convert SAM to BAM
if not os.path.exists(f'{work_dir}/need_second.bam'):
    cmd = f'{samtools} sort -O BAM --threads 8 {work_dir}/need_second.sam -o {work_dir}/need_second.bam'
    subprocess.run(cmd, shell=True, check=True)

# Index the BAM file
cmd = f'{samtools} index {work_dir}/need_second.bam'
subprocess.run(cmd, shell=True, check=True)
