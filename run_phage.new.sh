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
fastp -i $input_1 -I $input_2 -o $work_dir/${prefix}_1_filter.fq.gz -O $work_dir/${prefix}_2_filter.fq.gz -w 10

##prepared for seqGraph
bwa index $edge_fa
bwa mem -t 24 $edge_fa $work_dir/${prefix}_1_filter.fq.gz $work_dir/${prefix}_2_filter.fq.gz | samtools view -buSF 0x80c - > $work_dir/${prefix}_reads_pe_primary.bam
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
#python $find_phage_gene_matches -f $edge_fa -n 32 -o $work_dir -p /home/ruohawang2/12.Phage_assem/phage_protein_anno/SCAPP/scapp/protein_db

## node scoring
#python $phage_scoring $edge_fa $node_score True 32

## SeqGraph
samtools faidx $edge_fa
$rDistance -b $bam --out $work_dir/$prefix.txt  --depth $depth -t 0 -c -1
$extract_ref $fq1 $fq2 $phagedb $work_dir/${prefix}_ref_index2.txt 0.7 0.85 48 > $work_dir/${prefix}_ref_index.txt
echo "rDistance done"

python $get_ref_by_index $phagedb $phagedb.fai $work_dir/${prefix}_ref_index.txt $work_dir/${prefix}_ref_index.fa
makeblastdb -in $work_dir/${prefix}_ref_index.fa -dbtype nucl -parse_seqids -out $work_dir/${prefix}_ref_index

## Step3 recon and second iteration
mkdir -p $work_dir/second_match
cd $work_dir/second_match
python $recon $work_dir/$prefix.txt $prefix $need_second $depth
#deal the each  ${name}_{refname}.txt
for i in `ls $work_dir/second_match/*.second`; do
	fullname=$i
  	second=$(echo $fullname | sed 's/\.[^.]*$//');
  	$bin_m/matching -g $fullname -r ${second}_result.txt -c ${second}_cycle.txt -i 10 -v 1 --model 1 --ignore_copy;
  	python $make_fa_from_path $edge_fa ${second}_result.txt ${second}_unfiltered.fasta 1;
  	blastn -query ${second}_unfiltered.fasta -out ${second}_ass_index.blast -db $work_dir/${prefix}_ref_index -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" ;
  	python $plost ${second}_ass_index.blast $cycle_txt ${second}_tmp.txt 0 > ${second}_final.txt
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

## delete the unused files
#mkdir $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_final.fa $work_dir/${prefix}_final_result/
#mv $work_dir/filtered_cycle_res.txt $work_dir/${prefix}_final_result/
#mv $work_dir/need_second_match.txt $work_dir/${prefix}_final_result/
#mv $work_dir/${prefix}_assem/contigs.fasta $work_dir/${prefix}_final_result/

#mv $work_dir/${prefix}_final_result /home/ruohawang2/12.Phage_assem/meta_result/ERP003612/

rm $fq1
rm $fq2
