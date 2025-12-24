import re
import sys
from pathlib import Path

def parse_input_file(file_path, ignore_len):
    """Parse input file and extract valid lines."""
    res = set()
    
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if "loop" in line or "iter" in line:
                continue
            
            # Calculate line length if ignore_len is 0
            if ignore_len == 0:
                line_len = sum(
                    int(v.split('_')[3]) 
                    for v in re.split(r'[+-]', line) 
                    if v.strip()
                )
                if line_len < 10000:
                    continue
            
            # Clean the line
            cleaned_line = line
            for word in ["cycle", "score", "self", "gene", "ref"]:
                cleaned_line = cleaned_line.replace(word, "")
            
            res.add(cleaned_line.strip())
    
    return res

def load_gene_hits(gene_hit_file, min_count=5):
    """Load gene hits that meet minimum count threshold."""
    gene_hits = set()
    
    with open(gene_hit_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and int(parts[1]) >= min_count:
                gene_hits.add(parts[0])
    
    return gene_hits

def load_score_hits(score_file, min_score=0.7):
    """Load score hits that meet minimum score threshold."""
    score_hits = set()
    
    with open(score_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and float(parts[1]) >= min_score:
                score_hits.add(parts[0])
    
    return score_hits

#def filter_results(res, gene_hits, score_hits):
#    """Filter results based on gene and score hits."""
#    filtered_results = []
#    
#    for item in res:
#        # Extract contigs without direction indicators
#        contigs = [contig.rstrip('+-') for contig in re.split(r'(?=[+-])', item) if contig.strip()]
#        
#        # Single contig case: must be in gene_hits or score_hits
#        if len(contigs) <= 1:
#            if contigs and (contigs[0] in gene_hits or contigs[0] in score_hits):
#                filtered_results.append(item)
#        else:
#            # Multiple contigs: always include
#            filtered_results.append(item)
#    
#    return filtered_results

def filter_results(res, gene_hits, score_hits):
    """过滤结果并在每个contig（保留其正负号）后添加tab进行分割。"""
    filtered_results = []

    for item in res:
        # 使用findall提取每个以正负号结尾的片段
        contig_list = re.findall(r'.+?[+-]', item)
        # 用于判断：去掉各个片段末尾的正负号
        contigs_for_check = [contig.rstrip('+-') for contig in contig_list]

        # 单个contig情况：只有一个的情况，需要该contig在gene_hits或score_hits中
        if len(contigs_for_check) <= 1:
            if contigs_for_check and (contigs_for_check[0] in gene_hits or contigs_for_check[0] in score_hits):
                # 用\t拼接时，tab会在每个contig（其正负号之后）后面进行分割
                filtered_results.append("\t".join(contig_list))
        else:
            # 多个contigs情况始终保留
            filtered_results.append("\t".join(contig_list))

    return filtered_results
def main():
    if len(sys.argv) != 6:
        print("Usage: python script.py <input_file> <ignore_len> <gene_hit_file> <score_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    ignore_len = int(sys.argv[2])
    gene_hit_file = sys.argv[3]
    score_file = sys.argv[4]
    output_file = sys.argv[5]
    
    # Process input data
    res = parse_input_file(input_file, ignore_len)
    gene_hits = load_gene_hits(gene_hit_file)
    score_hits = load_score_hits(score_file)
    
    # Filter and write results
    filtered_results = filter_results(res, gene_hits, score_hits)
    
    with open(output_file, 'w') as f:
        for item in filtered_results:
            f.write(item + '\n')

if __name__ == "__main__":
    main()
