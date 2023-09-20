#!/bin/bash
assembly_fasta=$1
assembly_bam=$2
prefix=$3
work_dir=$4
threads=$5

##tools variable path
PYTHON=python
BWA=bwa
SAMTOOLS=~/apps/usr/local/bin/samtools
FASTP=~/apps/usr/local/bin/fastp
PALACE=/scratch/project/cs_shuaicli/pgz/phage/PALACE
bin_m="$PALACE/seqGraph_phage/build"
scripts="$PALACE/seqGraph_phage/scripts"
SPADES=/scratch/project/cs_shuaicli/pgz/phage/spades/bin/spades.py

make_fa_from_path="$scripts/make_fa_from_path.py"
make_fa_from_result="$scripts/make_fa_from_result.py"
make_fa_from_blast="$scripts/make_fa_from_blast.py"
make_final_fa="$scripts/make_final_fa.py"
filter_cycle="$scripts/filter_cycle.py"
get_ref_by_index="$scripts/get_ref_by_index.py"
plost="$scripts/plost.py"
recon="$scripts/recon.py"
make_fasta_from_fastg="$scripts/make_fasta_from_fastg.py"
find_phage_gene_matches="$scripts/find_phage_gene_matches.py"
phage_scoring="$scripts/phage_scoring.py"
rDistance="$bin_m/rDistance"
extract_ref="$bin_m/eref"
extract_fa_bwa="$scripts/extract_fa_bwa.py"
corrected_dup="$scripts/corrected_dup.py"
phagedb="/scratch/project/cs_shuaicli/pgz/phage/genome_db/checkv_reps.fna" #$4
protein_db="/scratch/project/cs_shuaicli/pgz/phage/protein_db"


hit_out="$work_dir/hit_seqs.out" #$7
node_score="$work_dir/node_scores.out" #$8

## Pre-assembly by metaspade
mkdir -p $work_dir
echo "Finish raw assembly"

# phage-related protein
if [ ! -s "$work_dir/hit_seqs.out" ]; then
	$PYTHON $find_phage_gene_matches -f $assembly_fasta -n $threads -o $work_dir -p ${protein_db}
fi
# deep_learning
if [ ! -s "$work_dir/node_scores.out" ]; then
	$PYTHON $phage_scoring $assembly_fasta $node_score True $threads
fi
# reference_based
# rm -r $work_dir


