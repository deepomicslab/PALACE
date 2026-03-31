import sys
import re
from collections import defaultdict

# ==========================================
# 1. Command line arguments
# ==========================================
fastg_fai_file = sys.argv[1]
original_graph_file = sys.argv[2]
output_file = sys.argv[3]
depth = int(float(sys.argv[4]))
f_th = sys.argv[5]
gene_file = sys.argv[6]
score_file = sys.argv[7]
blast_file = sys.argv[8]
blast_ratio = float(sys.argv[9])
fasta_fai_file = sys.argv[10]
hit_segs_file = sys.argv[11]
contig_path_file = sys.argv[12]
score_threshold = float(sys.argv[13])  # Score threshold parameter

# ==========================================
# 2. Global constants and data structures
# ==========================================
TO_REMOVE_SCORE_THRESHOLD = 0.2
RELEVATE_EDGE_LEN = 200
MIN_CYCLE_LEN = 1000
MIN_ALN_LEN = 1000
SAMPLE = 'SAMPLE'

num_to_full_name = {}
full_name_to_num = {}
fai_len = {}
gene_res = {}
scores = defaultdict(lambda: '0')
blast_segs = set()
relevate_blast_segs = set()  # Store 0-hop and 1-hop nodes
score_segs = set()
fastg_graph = {}
all_segs = {}
write_segs = set()
write_juncs = []
hit_segs = defaultdict(str)


# ==========================================
# 3. Helper parsing and filtering functions
# ==========================================
def get_edge_len(edge):
    """Extract length from node name"""
    return int(edge.split("_")[3])

def parse_fasta_index():
    """Parse FASTA index file to get sequence lengths and build name mapping"""
    with open(fasta_fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            seq_name = fields[0]
            seq_len = int(fields[1])
            fai_len[seq_name] = seq_len

            node_num = seq_name.split("_")[1]
            num_to_full_name[node_num] = seq_name
            full_name_to_num[seq_name] = node_num

def parse_blast_results():
    """Parse BLAST results to find nodes matching alignment ratio or length"""
    prev_seg = ""
    prev_ref = ""
    prev_len = 0

    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            query, ref, identity, aln_len = fields[0], fields[1], float(fields[2]), int(fields[3])

            if (prev_seg != query and prev_seg != "") or (prev_ref != ref and prev_ref != ""):
                seg_len = fai_len[prev_seg]
            # If alignment ratio meets threshold or length > 2000, add to blast_segs
                if prev_len / seg_len > blast_ratio or prev_len > 2000:
                    blast_segs.add(prev_seg)
                prev_seg = query
                prev_ref = ref
                prev_len = aln_len if identity > blast_ratio * 100 else 0
            else:
                if identity > blast_ratio * 100:
                    prev_len += aln_len
                prev_seg = query
                prev_ref = ref

    if prev_seg and prev_seg in fai_len:
        seg_len = fai_len[prev_seg]
        if prev_len / seg_len > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)

def parse_gene_and_score_files():
    """Parse gene annotation and score files"""
    if len(sys.argv) > 5:
        with open(gene_file, 'r') as f:
            for line in f:
                contig_name = line.split("\t")[0]
                gene_res[contig_name] = '1'

        with open(score_file, 'r') as f:
            for line in f:
                fields = line.strip().split("\t")
                contig_name, score_str = fields[0], fields[1]
                if 'e' in score_str.lower():
                    score_value = '0.0'
                else:
                    score_value = f"{float(score_str):.3f}"
                scores[contig_name] = score_value
                
                # Use passed score_threshold parameter instead of hardcoded 0.7
                if float(score_value) > score_threshold:
                    score_segs.add(contig_name)

def parse_fastg_index():
    """Parse FASTG index file"""
    with open(fastg_fai_file, 'r') as f:
        for line in f:
            fields = line.split("\t")
            parts = re.split(":|,|;", fields[0])
            fastg_graph[parts[0]] = parts[1:]

def filter_paths(support_segs):
    """Evaluate paths: if > half (or > 2000bp) consists of core seed, recover all nodes in path"""
    full_name_results = set()
    with open(contig_path_file, 'r') as f:
        for line in f:
            line = line.strip().replace(";", "")
            if line.startswith("NODE"):
                continue

            nums = line.split(",")
            full_names = []
            full_len = 0
            add_len = 0

            for num in nums:
                full_name = num_to_full_name[num[:-1]]
                full_names.append(full_name)
                e_len = get_edge_len(full_name)
                full_len += e_len
                if full_name in support_segs:
                    add_len += e_len

            if add_len > 0 and (add_len / full_len >= 0.5 or add_len > 2000):
                full_name_results.update(full_names)

    return full_name_results

def should_include_segment(seg_name):
    return (seg_name in blast_segs or
            seg_name in gene_res or
            float(scores.get(seg_name, '0')) > score_threshold)

def update_hit_segs(seg_name):
    hit_info = []
    if seg_name in blast_segs:
        hit_info.append("ref+")
        
    if float(scores.get(seg_name, '0')) > score_threshold:
        hit_info.append("score+")
        
    if seg_name in gene_res:
        hit_info.append("gene+")

    if hit_info:
        hit_segs[seg_name] = "".join(hit_info)
        relevate_blast_segs.add(seg_name)

def process_segment(seg_name, line):
    fields = line.strip().split()

    cleaned_fields = [fields[0], fields[1]]

    for field in fields[2:]:
        if 'e' in field.lower():
            try:
                val = float(field)
                if val.is_integer():
                    cleaned_fields.append(str(int(val)))
                else:
                    cleaned_fields.append(f"{val:.3f}".rstrip('0').rstrip('.'))
            except ValueError:
                cleaned_fields.append(field)
        else:
            cleaned_fields.append(field)

    cleaned_line = " ".join(cleaned_fields)

    is_blast = '1' if seg_name in blast_segs else '0'
    gene_val = gene_res.get(seg_name, '0')
    score_val = scores.get(seg_name, '0.000')

    return f"{cleaned_line} {gene_val} {score_val} {is_blast}\n"


# ==========================================
def main():
    parse_fasta_index()
    parse_blast_results()
    parse_gene_and_score_files()
    parse_fastg_index()

    with open(original_graph_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] == "SEG":
            seg_name = fields[1]
            all_segs[seg_name] = line
            update_hit_segs(seg_name)

            if should_include_segment(seg_name):
                write_segs.add(process_segment(seg_name, line))

    core_seeds = set(relevate_blast_segs) 
    hop1_segs = set()                     

    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] != "SEG":
            left_seg, right_seg = fields[1], fields[3]

            if left_seg == right_seg or left_seg in core_seeds or right_seg in core_seeds:
                write_juncs.append(line)
                write_segs.add(process_segment(left_seg, all_segs[left_seg]))
                write_segs.add(process_segment(right_seg, all_segs[right_seg]))
                hop1_segs.add(left_seg)
                hop1_segs.add(right_seg)

    relevate_blast_segs.update(hop1_segs)

    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] != "SEG":
            left_seg, right_seg = fields[1], fields[3]

            if left_seg in relevate_blast_segs or right_seg in relevate_blast_segs:
                write_juncs.append(line)
                write_segs.add(process_segment(left_seg, all_segs[left_seg]))
                write_segs.add(process_segment(right_seg, all_segs[right_seg]))

    support_segs = blast_segs | set(gene_res.keys()) | score_segs
    path_segs = filter_paths(support_segs)
    written_segs = {item.split(" ")[1] for item in write_segs}

    # ==========================================
    with open(output_file, 'w') as out:
        for seg_line in write_segs:
            out.write(seg_line)

        for seg in path_segs:
            if seg not in written_segs:
                out.write(f"{all_segs[seg].strip()} 0 1.0 0\n")

        seen_juncs = set()
        for junc in write_juncs:
            if junc not in seen_juncs:
                out.write(junc)
                seen_juncs.add(junc)

    with open(hit_segs_file, 'w') as out:
        for seg_name, hit_info in hit_segs.items():
            if hit_info:
                out.write(f"{SAMPLE}\t{seg_name}\t{hit_info}\n")


if __name__ == "__main__":
    main()
