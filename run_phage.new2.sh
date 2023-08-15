#!/bin/bash
input_1=$1
input_2=$2
prefix=$3
work_dir=$4
edge_fa=$5
cycle_txt=$6
need_second=$7

##tools variable path
bin_m="/home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build"
scripts="/home/ruohawang2/12.Phage_assem/pipeline/seqGraph/scripts"

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

## Step1
##raw assembly by metaspade
bam="$work_dir/${prefix}_reads_pe_primary.sort.bam"  #$1
fq1="$work_dir/${prefix}_1_filter.fq" #$2
fq2="$work_dir/${prefix}_2_filter.fq" #$3
phagedb="/home/ruohawang2/12.Phage_assem/pipeline/genome_db/checkv_reps.fna" #$4
assembly_fasta="$work_dir/${prefix}_assem/assembly_graph.fasta" #$5
assembly_fastg="$work_dir/${prefix}_assem/assembly_graph.fastg" #$6

hit_out="$work_dir/hit_seqs.out" #$7
node_score="$work_dir/node_scores.out" #$8
out_prefix=$prefix #$9

## protein annotation
#python $find_phage_gene_matches -f $edge_fa -n 32 -o $work_dir -p /home/ruohawang2/12.Phage_assem/phage_protein_anno/SCAPP/scapp/protein_db

## node scoring
#python $phage_scoring $edge_fa $node_score True 32

## SeqGraph
#python $filter_cycle $work_dir/cycle_res.txt > $work_dir/filtered_cycle_res.txt
if [ -f "$cycle_txt" ]; then
	cat $cycle_txt > $work_dir/${prefix}_final.txt
fi

if [ "$(ls -A $work_dir/second_match/*_final.txt 2>/dev/null)" ]; then
	cat $work_dir/second_match/*_final.txt >> $work_dir/${prefix}_final.txt
fi

if [ "$(ls -A $work_dir/second_match/*cycle.txt 2>/dev/null)" ]; then
	cat $work_dir/second_match/*cycle.txt | grep -v 'iter' >> $work_dir/${prefix}_final.txt
fi

echo "second match done\n"

python $make_final_fa $edge_fa $work_dir/${prefix}_final.txt $work_dir/${prefix}_final.fa $prefix

## delete the unused files
#mkdir $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.fa $work_dir/${prefix}_final_result/
#mv $work_dir/filtered_cycle_res.txt $work_dir/${prefix}_final_result/
#mv $work_dir/need_second_match.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_assem/contigs.fasta $work_dir/${prefix}_final_result/

#mv $work_dir/${prefix}_final_result /home/ruohawang2/12.Phage_assem/meta_result/ERP003612/


