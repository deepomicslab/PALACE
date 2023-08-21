import subprocess
import re
import sys
# Define the files and the working directory
file_a = sys.argv[1]
file_b = sys.argv[2]
fq1 = sys.argv[3]
fq2 = sys.argv[4]
work_dir = sys.argv[5]

# Read file A and create the string of regions to extract
line_all = []
with open(file_a, 'r') as f:
    for line in f:
        line_arr = re.split('[-+\t]', line.strip())
        line_all.append(' '.join(line_arr[:-1]))
line_str = ' '.join(line_all)

# Extract the subfasta
cmd = f'samtools faidx {file_b} {line_str} > {work_dir}/need_second.fasta'
subprocess.run(cmd, shell=True, check=True)

# Index the subfasta
cmd = f'bwa index {work_dir}/need_second.fasta'
subprocess.run(cmd, shell=True, check=True)

# Align the fastq files to the subfasta
cmd = f'bwa mem 8 {work_dir}/need_second.fasta {fq1} {fq2} > {work_dir}/need_second.sam'
subprocess.run(cmd, shell=True, check=True)

# Convert SAM to BAM
cmd = f'samtools view -bS {work_dir}/need_second.sam > {work_dir}/need_second.bam'
subprocess.run(cmd, shell=True, check=True)

# Index the BAM file
cmd = f'samtools index {work_dir}/need_second.bam'
subprocess.run(cmd, shell=True, check=True)