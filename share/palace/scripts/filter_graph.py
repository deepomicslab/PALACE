import sys
import re
from collections import defaultdict

# Parse command line arguments
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

# Constants
TO_REMOVE_SCORE_THRESHOLD = 0.2
RELEVATE_EDGE_LEN = 200
MIN_CYCLE_LEN = 1000
MIN_ALN_LEN = 1000
SAMPLE = 'SAMPLE'

# Initialize data structures
num_to_full_name = {}
full_name_to_num = {}
fai_len = {}
gene_res = {}
scores = defaultdict(lambda: '0')
blast_segs = set()
relevate_blast_segs = set()
score_segs = set()
fastg_graph = {}
all_segs = {}
write_segs = set()
write_juncs = []
hit_segs = defaultdict(str)

# Helper functions
def get_edge_len(edge):
    """Extract edge length from edge name"""
    return int(edge.split("_")[3])

def parse_fasta_index():
    """Parse FASTA index file to get sequence lengths and name mappings"""
    with open(fasta_fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            seq_name = fields[0]
            seq_len = int(fields[1])
            fai_len[seq_name] = seq_len
            
            # Extract node number from sequence name
            node_num = seq_name.split("_")[1]
            num_to_full_name[node_num] = seq_name
            full_name_to_num[seq_name] = node_num

def parse_blast_results():
    """Parse BLAST results to identify segments with significant hits"""
    prev_seg = ""
    prev_ref = ""
    prev_len = 0
    
    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            query, ref, identity, aln_len = fields[0], fields[1], float(fields[2]), int(fields[3])
            
            # Check if we've moved to a new query or reference
            if (prev_seg != query and prev_seg != "") or (prev_ref != ref and prev_ref != ""):
                # Check if previous segment meets criteria
                seg_len = fai_len[prev_seg]
                if prev_len / seg_len > blast_ratio or prev_len > 2000:
                    blast_segs.add(prev_seg)
                prev_seg = query
                prev_ref = ref
                prev_len = aln_len if identity > blast_ratio * 100 else 0
            else:
                # Accumulate alignment length for high-identity hits
                if identity > blast_ratio * 100:
                    prev_len += aln_len
                prev_seg = query
                prev_ref = ref
    
    # Process the last segment
    if prev_seg and prev_seg in fai_len:
        seg_len = fai_len[prev_seg]
        if prev_len / seg_len > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)

def parse_gene_and_score_files():
    """Parse gene and score files"""
    # Parse gene file
    if len(sys.argv) > 5:
        with open(gene_file, 'r') as f:
            for line in f:
                contig_name = line.split("\t")[0]
                gene_res[contig_name] = '1'
        
        # Parse score file
        with open(score_file, 'r') as f:
            for line in f:
                fields = line.strip().split("\t")
                contig_name, score_value = fields[0], fields[1]
                scores[contig_name] = score_value
                if float(score_value) > 0.7:
                    score_segs.add(contig_name)

def parse_fastg_index():
    """Parse FASTG index file"""
    with open(fastg_fai_file, 'r') as f:
        for line in f:
            fields = line.split("\t")
            parts = re.split(":|,|;", fields[0])
            fastg_graph[parts[0]] = parts[1:]

def filter_paths(support_segs):
    """Filter contig paths based on support segments"""
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
            
            # Check if path meets criteria
            #if add_len > 0 and (add_len / full_len >= 0.5 or add_len > 2000):
            full_name_results.update(full_names)
    
    return full_name_results

def process_segment(seg_name, line):
    """Process a segment and return formatted output"""
    is_blast = '1' if seg_name in blast_segs else '0'
    gene_val = gene_res.get(seg_name, '0')
    score_val = scores.get(seg_name, '0')
    return f"{line.strip()} {gene_val} {score_val} {is_blast}\n"

def should_include_segment(seg_name):
    """Check if segment should be included based on various criteria"""
    return (seg_name in blast_segs or 
            seg_name in gene_res or 
            float(scores.get(seg_name, '0')) > 0.7 or
            seg_name in relevate_blast_segs)

def update_hit_segs(seg_name):
    """Update hit segments tracking"""
    hit_info = []
    if seg_name in blast_segs:
        hit_info.append("ref+")
    if float(scores.get(seg_name, '0')) > 0.7:
        hit_info.append("score+")
    if seg_name in gene_res:
        hit_info.append("gene+")
    
    if hit_info:
        hit_segs[seg_name] = "".join(hit_info)
        relevate_blast_segs.add(seg_name)

# Main processing
def main():
    # Parse all input files
    parse_fasta_index()
    parse_blast_results()
    parse_gene_and_score_files()
    parse_fastg_index()
    
    # Read original graph
    with open(original_graph_file, 'r') as f:
        lines = f.readlines()
    
    # First pass: process segments and identify relevant segments
    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] == "SEG":
            seg_name = fields[1]
            all_segs[seg_name] = line
            update_hit_segs(seg_name)
            
            if should_include_segment(seg_name):
                write_segs.add(process_segment(seg_name, line))
    
    # Process junctions (first criteria)
    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] != "SEG":
            left_seg, right_seg = fields[1], fields[3]
            left_score = float(scores.get(left_seg, '0'))
            right_score = float(scores.get(right_seg, '0'))
            
            # Check various conditions for including junction
            if (left_seg == right_seg or
                should_include_segment(left_seg) or 
                should_include_segment(right_seg)):
                
                write_juncs.append(line)
                write_segs.add(process_segment(left_seg, all_segs[left_seg]))
                write_segs.add(process_segment(right_seg, all_segs[right_seg]))
                relevate_blast_segs.add(left_seg)
                relevate_blast_segs.add(right_seg)
    
    # Second pass: include junctions connected to relevant segments
    for line in lines:
        fields = line.rstrip().split(" ")
        if fields[0] != "SEG":
            left_seg, right_seg = fields[1], fields[3]
            if left_seg in relevate_blast_segs or right_seg in relevate_blast_segs:
                write_juncs.append(line)
                write_segs.add(process_segment(left_seg, all_segs[left_seg]))
                write_segs.add(process_segment(right_seg, all_segs[right_seg]))
    
    # Filter paths and add missing segments
    support_segs = blast_segs | set(gene_res.keys()) | score_segs
    path_segs = filter_paths(support_segs)
    
    # Extract already written segments
    written_segs = {item.split(" ")[1] for item in write_segs}
    
    # Write output
    with open(output_file, 'w') as out:
        # Write segments
        for seg_line in write_segs:
            out.write(seg_line)
        
        # Write path segments not already written
        for seg in path_segs:
            if seg not in written_segs:
                out.write(f"{all_segs[seg].strip()} 0 1.0 0\n")
        
        # Write unique junctions
        seen_juncs = set()
        for junc in write_juncs:
            if junc not in seen_juncs:
                out.write(junc)
                seen_juncs.add(junc)
    
    # Write hit segments file
    with open(hit_segs_file, 'w') as out:
        for seg_name, hit_info in hit_segs.items():
            if hit_info:
                out.write(f"{SAMPLE}\t{seg_name}\t{hit_info}\n")

if __name__ == "__main__":
    main()
