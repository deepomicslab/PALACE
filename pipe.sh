#!/bin/bash

###############################parse config file###############################
config_file=$1
while IFS='=' read -r key value
do
    # Skip lines starting with #
    [[ "$key" =~ ^\#.* ]] || [[ -z "$key" ]] && continue
    key=$(echo $key | tr '.' '_')
    echo $key"="$value
    eval "${key}='${value}'"
done < "$config_file"
#print time
print_time() {
    echo $(date +"%Y-%m-%d %H:%M:%S")
}
#creat directory
create_dir() {
    dir_name=$1
    if [ ! -d "$dir_name" ]; then
    mkdir -p "$dir_name"
    echo "Directory $dir_name has been created."
    else
    echo "Directory $dir_name already exists."
    fi
}
create_dir $out_dir
###############################SCRIPTS###############################
MATCHING=$PALACE/seqGraph_phage/build/matching
RDistance=$PALACE/seqGraph_phage/build/rDistance
EXTRACT_REF=$PALACE/seqGraph_phage/build/eref

SCRIPTS=$PALACE/seqGraph_phage/scripts
MAKE_FA_FROM_PATH="$SCRIPTS/make_fa_from_path.py"
FILTER_RESULT="$SCRIPTS/filter_result.py"
MAKE_FINAL_FA="$SCRIPTS/make_final_fa.py"
FILTER_CYCLE="$SCRIPTS/filter_cycle.py"
GET_REF_BY_INDEX="$SCRIPTS/get_ref_by_index.py"
FILTER_BY_BLAST="$SCRIPTS/filter_by_blast.py"
EXTRACT_BY_REF="$SCRIPTS/extract_by_ref.py"
SPLIT_FASTG="$SCRIPTS/split_fastg.py"
FIND_PHAGE_GENE_MATCHES="$SCRIPTS/find_phage_gene_matches.py"
PHAGE_SCORING="$SCRIPTS/phage_scoring.py"
CORRECTED_DUP="$SCRIPTS/corrected_dup.py"
FILTER_GRAPH="$SCRIPTS/filter_graph.py"

###############################intermediate files###############################
filter_fastq1="$out_dir/01-qc/${prefix}_1_filter.fastq"
filter_fastq2="$out_dir/01-qc/${prefix}_2_filter.fastq"
first_bam="$out_dir/02-assembly/${prefix}_reads_pe_primary.sort.bam"
assembly_fasta="$out_dir/02-assembly/assembly_graph.fasta"
assembly_fastg="$out_dir/02-assembly/assembly_graph.fastg"
hit_out="$out_dir/03-search/hit_seqs.out"
node_score="$out_dir/03-search/node_scores.out"
phage_refs="$out_dir/03-search/phage_refs.fasta"


###############################Step 1, fastqc###############################
create_dir $out_dir/01-qc
echo "$(print_time) Step 1, fastq QC..."
if [[ ! -s "$filter_fastq1" && ! -s "$filter_fastq2" ]]; then
        $FASTP -i $fastq1 -I $fastq2 -o $filter_fastq1 -O $filter_fastq2 -w $threads
fi
echo "$(print_time) Finished Step 1"

###############################Step 2, Raw assembly and align###############################
create_dir $out_dir/02-assembly
echo "$(print_time) Step 2, assembly"
if [ ! -s "$out_dir/02-assembly/contigs.fasta" ]; then
    $SPADES --meta -o $out_dir/02-assembly -1 $filter_fastq1 -2 $filter_fastq2 -t $threads --phred-offset 33
fi
$PYTHON $SPLIT_FASTG -g $assembly_fastg -o $assembly_fasta
echo "$(print_time) Step 2, assembly done"
echo "$(print_time) Step 2, start align..."
$SAMTOOLS faidx $assembly_fasta
$SAMTOOLS faidx $assembly_fastg
if [ ! -s "$first_bam.bai" ]; then
    $BWA index $assembly_fasta
    $BWA mem -t $threads $assembly_fasta $filter_fastq1 $filter_fastq2 | $SAMTOOLS view -F 0x0800 -buS - > $out_dir/02-assembly/${prefix}_tmp.bam
    $SAMTOOLS sort -@ $threads $out_dir/02-assembly/${prefix}_tmp.bam -O BAM -o $first_bam
    $SAMTOOLS index $first_bam
fi
echo "$(print_time) Step 2, align done"
first_depth=`$SAMTOOLS depth $first_bam | awk '{sum+=$3} END { print sum/NR}'` #$10
echo "$(print_time) Finished Step 2"

###############################Step 3, search ref from database, search protein and predict phage contigs###############################
create_dir $out_dir/03-search
echo "$(print_time) Step 3, searching protein..."
if [ ! -s "$hit_out" ]; then
    $PYTHON $FIND_PHAGE_GENE_MATCHES -f $assembly_fasta -n $threads -o $out_dir/03-search -p ${protein_db}
fi
echo "$(print_time) Step 3, search protein done"
echo "$(print_time) Step 3, searching contigs..."
# deep_learning
if [ ! -s "$node_score" ]; then
        echo "$PYTHON $PHAGE_SCORING $assembly_fasta $node_score True $threads $gcn_model"
    $PYTHON $PHAGE_SCORING $assembly_fasta $node_score True $threads $gcn_model
fi
echo "$(print_time) Step 3, search contigs done"
echo "$(print_time) Step 3, searching reference..."
if [ ! -f "$out_dir/03-search/${prefix}_ref_names.txt" ]; then
    $EXTRACT_REF $filter_fastq1 $filter_fastq2 $phagedb $out_dir/03-search/${prefix}_tmp.txt 0.9 0.85 48 > $out_dir/03-search/${prefix}_ref_names.txt
fi
if [ ! -f "$phage_refs" ]; then
    $PYTHON $GET_REF_BY_INDEX $phagedb $phagedb.fai $out_dir/03-search/${prefix}_ref_names.txt $phage_refs
fi
echo "$(print_time) Step 3, searching done"
echo "$(print_time) Finished Step 3"

###############################Step 4, constract graph and matching cycles###############################
create_dir $out_dir/04-match
echo "$(print_time) Step 4, creating graph..."
if [ ! -s "$assembly_fasta.blast" ]; then
    $NCBI_BIN/makeblastdb -in $phage_refs -dbtype nucl -parse_seqids -out $phage_refs
    blastn -query $assembly_fasta -out $assembly_fasta.blast -db $phage_refs -outfmt 6
fi
if [ ! -f "$out_dir/04-match/${prefix}_graph.txt" ]; then
    $RDistance -b $first_bam --out $out_dir/04-match/${prefix}_graph.txt --depth $first_depth -t 0 -c -1
fi
if [ ! -f "$out_dir/04-match/${prefix}_filtered_graph.txt" ]; then
## Collect the filtered contigs, meeting one requirement, and the connected contigs
    $PYTHON $FILTER_GRAPH $assembly_fastg.fai $out_dir/04-match/${prefix}_graph.txt $out_dir/04-match/${prefix}_filtered_graph.txt $first_depth 0 $hit_out $node_score $assembly_fasta.blast  0.7 $assembly_fasta.fai
fi
echo "$(print_time) Step 4, graph created"
echo "$(print_time) Step 4, start matching"
## Maximum matching, outputing linear and circular sequences
$MATCHING -g $out_dir/04-match/${prefix}_filtered_graph.txt -r $out_dir/04-match/${prefix}_linear.txt -c $out_dir/04-match/${prefix}_cycle.txt -i 10 -v 1 -s
cat $out_dir/04-match/${prefix}_linear.txt $out_dir/04-match/${prefix}_cycle.txt > $out_dir/04-match/${prefix}_all_result.txt
# filter
if [ ! -f "$out_dir/04-match/${prefix}_unfiltered.fasta" ];then
    $PYTHON $MAKE_FA_FROM_PATH $assembly_fasta $out_dir/04-match/${prefix}_all_result.txt $out_dir/04-match/${prefix}_unfiltered.fasta 0
fi
if [ ! -f "$out_dir/04-match/${prefix}_unfiltered_phagescore.txt" ];then
    $PYTHON $PHAGE_SCORING $out_dir/04-match/${prefix}_unfiltered.fasta $out_dir/04-match/${prefix}_unfiltered_phagescore.txt False $threads $gcn_model
fi
if [ ! -f "$out_dir/04-match/${prefix}_filtered.fasta" ];then
    $PYTHON $FILTER_RESULT $assembly_fasta $out_dir/04-match/${prefix}_all_result.txt $out_dir/04-match/${prefix}_filtered.fasta $assembly_fasta.blast 0.75 $hit_out $out_dir/04-match/${prefix}_unfiltered_phagescore.txt $out_dir/04-match/${prefix}_filtered_cycle.txt
fi
if [ ! -f "$out_dir/04-match/${prefix}_filtered.fasta.blast" ];then
    blastn -query $out_dir/04-match/${prefix}_filtered.fasta -out $out_dir/04-match/${prefix}_filtered.fasta.blast -db $phage_refs -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"
fi
echo "Finished Step 4"

###############################Step 5, furth assembly with references###############################
create_dir $out_dir/05-furth/second_match
echo "Starting Step 5, furth assemblying..."
## Meet two requirements
$PYTHON $FILTER_BY_BLAST $out_dir/04-match/${prefix}_filtered.fasta.blast $out_dir/04-match/${prefix}_filtered_cycle.txt $assembly_fasta.fai $out_dir/05-furth/need_second_match.txt 0 0.7 2000 > $out_dir/04-match/${prefix}_tmp_cycle.txt
contig_names=$(cut -f 2 $out_dir/05-furth/need_second_match.txt | tr '\n' ' ')
$SAMTOOLS faidx ${phage_refs} ${contig_names} > $out_dir/05-furth/need_second.fasta
#sub_depth=`$SAMTOOLS depth ${first_bam} ${contig_names} |  awk '{sum+=$3} END { print sum/NR}'` #$10
$PYTHON $EXTRACT_BY_REF $out_dir/04-match/${prefix}_graph.txt $out_dir/05-furth/second_match/$prefix $out_dir/05-furth/need_second_match.txt ${SAMTOOLS} 5 ${first_bam} 0.7
# $PYTHON $recon $assembly_fastg.fai $out_dir/${prefix}_ass_blast.txt $prefix $out_dir/need_second_match.txt $depth
echo "recon done, Start Second Matching\n"

# iteraction: deal the each  ${name}_{refname}.txt
for i in `ls $out_dir/05-furth/second_match/*.second`; do
    fullname=$i
    second=$(echo $fullname | sed 's/\.[^.]*$//');
    $MATCHING -g $fullname -r ${second}_linear.txt -c ${second}_cycle.txt -i 10 -v 1 -b;
    if [ $? -eq 124 ]
    then
        echo "The command was terminated because it ran for more than 30 minutes."
    fi
    echo $second" second match"
    $PYTHON $MAKE_FA_FROM_PATH $assembly_fasta ${second}_linear.txt ${second}_unfiltered.fasta 1;
    blastn -query ${second}_unfiltered.fasta -out ${second}_unfiltered.fasta.blast -db $phage_refs -outfmt "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" ;
    $PYTHON $FILTER_BY_BLAST ${second}_unfiltered.fasta.blast ${second}_cycle.txt $assembly_fasta.fai ${second}_tmp.txt 0 0.7 2000 ${second}> ${second}_all_result.txt
done
echo "$(print_time) Finished Step 5"

###############################Step 6, make final result###############################
create_dir $out_dir/final_result
$PYTHON $FILTER_CYCLE $out_dir/04-match/${prefix}_filtered_cycle.txt 1 > $out_dir/final_result/filtered_cycle_res_tmp.txt
if [ -f "$out_dir/final_result/filtered_cycle_res_tmp.txt" ]; then
    cat $out_dir/final_result/filtered_cycle_res_tmp.txt > $out_dir/final_result/${prefix}_final_tmp.txt
fi

if [ "$(ls -A $out_dir/05-furth/second_match/*_all_result.txt 2>/dev/null)" ]; then
    cat $out_dir/05-furth/second_match/*_all_result.txt >> $out_dir/final_result/${prefix}_final_tmp.txt
fi

if [ "$(ls -A $out_dir/05-furth/second_match/*cycle.txt 2>/dev/null)" ]; then
    cat $out_dir/05-furth/second_match/*cycle.txt | grep -v 'iter' >> $out_dir/final_result/${prefix}_final_tmp.txt
fi

$PYTHON $MAKE_FINAL_FA $assembly_fasta $out_dir/final_result/${prefix}_final_tmp.txt $out_dir/final_result/${prefix}_final_tmp.fasta $prefix
$PYTHON $CORRECTED_DUP $out_dir/final_result ${prefix} $out_dir/final_result/filtered_cycle_res_tmp.txt $out_dir/final_result/${prefix}_final_tmp.txt ${prefix}_final.txt ${prefix}_final.fasta $assembly_fasta ${prefix}_cycle.txt $first_bam
echo "all done"