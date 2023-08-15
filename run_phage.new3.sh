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
mkdir -p $work_dir
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
mkdir -p $work_dir/second_match
cd $work_dir/second_match
blastn -query $edge_fa -out $work_dir/second_match/second_blast.txt -db $work_dir/${prefix}_ref_index -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"
samtools index $bam
python $recon $work_dir/$prefix.txt $prefix $need_second -1 5 $bam $work_dir/second_match/second_blast.txt 0.7
#deal the each  ${name}_{refname}.txt
for i in `ls $work_dir/second_match/*.second`; do
	fullname=$i
  	second=$(echo $fullname | sed 's/\.[^.]*$//');
  	$bin_m/matching -g $fullname -r ${second}_result.txt -c ${second}_cycle.txt -i 10 -v 1 --model 1;
  	python $make_fa_from_path $edge_fa ${second}_result.txt ${second}_unfiltered.fasta 1;
  	blastn -query ${second}_unfiltered.fasta -out ${second}_ass_index.blast -db $work_dir/${prefix}_ref_index -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" ;
  	python $plost ${second}_ass_index.blast $cycle_txt ${second}_tmp.txt 0 0.7> ${second}_final.txt
done

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
