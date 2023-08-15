#!/bin/bash
input_1=$1
input_2=$2
prefix=$3
work_dir=$4

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
mkdir -p $work_dir
fastp -i $input_1 -I $input_2 -o $work_dir/${prefix}_1_filter.fq.gz -O $work_dir/${prefix}_2_filter.fq.gz -w 32
spades.py --meta -o $work_dir/${prefix}_assem -1 $work_dir/${prefix}_1_filter.fq.gz -2 $work_dir/${prefix}_2_filter.fq.gz -t 32 --phred-offset 33
python $make_fasta_from_fastg -g $work_dir/${prefix}_assem/assembly_graph.fastg -o $work_dir/${prefix}_assem/assembly_graph.fasta
echo "Finish raw assembly"

##prepared for seqGraph
bwa index $work_dir/${prefix}_assem/assembly_graph.fasta
bwa mem -t 32 $work_dir/${prefix}_assem/assembly_graph.fasta $work_dir/${prefix}_1_filter.fq.gz $work_dir/${prefix}_2_filter.fq.gz | samtools view -buS - > $work_dir/${prefix}_reads_pe.bam
samtools view -bF 0x0800 $work_dir/${prefix}_reads_pe.bam > $work_dir/${prefix}_reads_pe_primary.bam
samtools sort -@ 8 $work_dir/${prefix}_reads_pe_primary.bam > $work_dir/${prefix}_reads_pe_primary.sort.bam
echo "Finish BWA"


## Step2
gunzip $work_dir/${prefix}_1_filter.fq.gz
gunzip $work_dir/${prefix}_2_filter.fq.gz
bam="$work_dir/${prefix}_reads_pe_primary.sort.bam"  #$1
fq1="$work_dir/${prefix}_1_filter.fq" #$2
fq2="$work_dir/${prefix}_2_filter.fq" #$3
phagedb="/home/ruohawang2/12.Phage_assem/pipeline/genome_db/checkv_reps.fna" #$4
assembly_fasta="$work_dir/${prefix}_assem/assembly_graph.fasta" #$5
assembly_fastg="$work_dir/${prefix}_assem/assembly_graph.fastg" #$6

hit_out="$work_dir/hit_seqs.out" #$7
node_score="$work_dir/node_scores.out" #$8
out_prefix=$prefix #$9
depth=`samtools depth $bam  |  awk '{sum+=$3} END { print sum/NR}'` #$10

## protein annotation
python $find_phage_gene_matches -f $assembly_fasta -n 32 -o $work_dir -p /home/ruohawang2/12.Phage_assem/phage_protein_anno/SCAPP/scapp/protein_db

## node scoring
python $phage_scoring $assembly_fasta $node_score True 32

## SeqGraph
samtools faidx $assembly_fasta
samtools faidx $assembly_fastg
$rDistance -b $bam --out $work_dir/$prefix.txt  --depth $depth -t 0 -c -1
$extract_ref $fq1 $fq2 $phagedb $work_dir/${prefix}_ref_index2.txt 0.7 0.85 48 > $work_dir/${prefix}_ref_index.txt
echo "rDistance done"

python $get_ref_by_index $phagedb $phagedb.fai $work_dir/${prefix}_ref_index.txt $work_dir/${prefix}_ref_index.fa
makeblastdb -in $work_dir/${prefix}_ref_index.fa -dbtype nucl -parse_seqids -out $work_dir/${prefix}_ref_index
blastn -query $assembly_fasta -out $work_dir/${prefix}_ref_index.blast -db $work_dir/${prefix}_ref_index -outfmt 6

python $scripts/parse.py $assembly_fastg.fai $work_dir/$prefix.txt $work_dir/${prefix}_ass_blast.txt 40 0 $hit_out $node_score $work_dir/${prefix}_ref_index.blast 0.7
$bin_m/matching -g $work_dir/${prefix}_ass_blast.txt -r $work_dir/${prefix}_result.txt -c $work_dir/${prefix}_cycle.txt -i 10 -v 1 -s
cat $work_dir/${prefix}_result.txt $work_dir/${prefix}_cycle.txt > $work_dir/${prefix}_all_r.txt
echo "Finish matching"

python $make_fa_from_path $assembly_fasta $work_dir/${prefix}_all_r.txt $work_dir/${prefix}_unfiltered.fasta 0
python $phage_scoring $work_dir/${prefix}_unfiltered.fasta $work_dir/${prefix}_unfiltered_phagescore.txt False 32
python $make_fa_from_result $assembly_fasta $work_dir/${prefix}_all_r.txt $work_dir/${prefix}_blast_seg.fasta $work_dir/${prefix}_ref_index.blast 0.75 $hit_out $work_dir/${prefix}_unfiltered_phagescore.txt $work_dir/cycle_res.txt
blastn -query $work_dir/${prefix}_blast_seg.fasta -out $work_dir/${prefix}_ass_index.blast -db $work_dir/${prefix}_ref_index  -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"

##plots
python $plost $work_dir/${prefix}_ass_index.blast $work_dir/cycle_res.txt $work_dir/need_second_match.txt 0 > $work_dir/${prefix}_tmp_cycle.txt


## Step3 recon and second iteration
mkdir -p $work_dir/second_match
cd $work_dir/second_match
python $recon $assembly_fastg.fai $work_dir/${prefix}_ass_blast.txt $prefix $work_dir/need_second_match.txt $depth
echo "recon done, Start Second Matching\n"

#deal the each  ${name}_{refname}.txt
if [ -e "$work_dir/second_match/*.second" ]; then
	for i in `ls $work_dir/second_match/*.second`; do
  	  fullname=$i
  	  second=$(echo $fullname | sed 's/\.[^.]*$//');
  	  $bin_m/matching -g $fullname -r ${second}_result.txt -c ${second}_cycle.txt -i 10 -v 1 -b --model 1 --ignore_copy;
  	  python $make_fa_from_path $assembly_fasta ${second}_result.txt ${second}_unfiltered.fasta 1;
  	  blastn -query ${second}_unfiltered.fasta -out ${second}_ass_index.blast -db $work_dir/${prefix}_ref_index -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" ;
  	  python $plost ${second}_ass_index.blast $work_dir/cycle_res.txt ${second}_tmp.txt 0 > ${second}_final.txt
	done
fi

python $filter_cycle $work_dir/cycle_res.txt > $work_dir/filtered_cycle_res.txt

if [ -e "$work_dir/second_match/*_final.txt" ]; then
	cat $work_dir/filtered_cycle_res.txt $work_dir/second_match/*_final.txt > $work_dir/${prefix}_final.txt
else
	cat $work_dir/filtered_cycle_res.txt > $work_dir/${prefix}_final.txt
fi

echo "second match done\n"


python $make_final_fa $assembly_fasta $work_dir/${prefix}_final.txt $work_dir/${prefix}_final.fa $prefix

## delete the unused files
#mkdir $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.fa $work_dir/${prefix}_final_result/
#mv $work_dir/filtered_cycle_res.txt $work_dir/${prefix}_final_result/
#mv $work_dir/need_second_match.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_assem/contigs.fasta $work_dir/${prefix}_final_result/

#mv $work_dir/${prefix}_final_result /home/ruohawang2/12.Phage_assem/meta_result/ERP003612/
#rm -r $work_dir


