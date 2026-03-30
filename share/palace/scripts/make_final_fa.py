#!/usr/bin/env python3
import sys
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_graph(graph_file):

    adjacency = {}

    with open(graph_file, 'r') as f:
        for line in f:
            if line.startswith('JUNC'):
                parts = line.strip().split()
                if len(parts) >= 5:
                    edge1, orient1, edge2, orient2 = parts[1:5]

                    src_node = f"{edge1}{orient1}"
                    dst_node = f"{edge2}{orient2}"

                    # Add direct connection
                    if src_node not in adjacency:
                        adjacency[src_node] = set()
                    adjacency[src_node].add(dst_node)

                    # Add conjugate connection (reverse complement)
                    conj_src = f"{edge2}{'+' if orient2 == '-' else '-'}"
                    conj_dst = f"{edge1}{'+' if orient1 == '-' else '-'}"

                    if conj_src not in adjacency:
                        adjacency[conj_src] = set()
                    adjacency[conj_src].add(conj_dst)

    return adjacency

def get_length_from_name(node_name):
    """Extract length from node name like EDGE_4484859_length_280_cov_+"""
    match = re.search(r'length_(\d+)', node_name)
    if match:
        return int(match.group(1))
    return float('inf') 

def is_circular_path_fuzzy(path, adjacency, trim_threshold, min_cycle_length):
    """
    Check if path forms a cycle (Fuzzy version).
    Constraint 1: Total length of trimmed contigs on both ends <= trim_threshold.
    Constraint 2: Total unique contig length in remaining cycle >= min_cycle_length.
    """
    if not path:
        return False, []

    # Pre-fetch lengths of all nodes in path
    lengths = [get_length_from_name(node) for node in path]
    valid_cycles = []

    # Iterate over all possible retention intervals [i, j]
    for i in range(len(path)):
        for j in range(i, len(path)):
            # Calculate total length of trimmed nodes on both ends
            trimmed_length = sum(lengths[:i]) + sum(lengths[j+1:])
            
            # Constraint 1: if trimmed length exceeds threshold, skip
            if trimmed_length > trim_threshold:
                continue
                
            first_node = path[i]
            last_node = path[j]
            
            # Check if stripped ends can connect in graph to form cycle
            if last_node in adjacency and first_node in adjacency[last_node]:
                subpath = path[i:j+1]
                
                # Get unique base contigs in cycle core (strip +/-)
                unique_base_edges = {node.rstrip('+-') for node in subpath}
                
                # Calculate total physical length after deduplication
                cycle_physical_length = sum(get_length_from_name(edge) for edge in unique_base_edges)
                
                # Constraint 2: if real cycle length < minimum, skip
                if cycle_physical_length >= min_cycle_length:
                    valid_cycles.append((trimmed_length, subpath))

    # If valid cycles found
    if valid_cycles:
        # Sort by trimmed length (ascending) to get the most complete cycle
        valid_cycles.sort(key=lambda x: x[0])
        return True, valid_cycles[0][1]

    return False, []

def make_fa_from_paths(fain, cycle_paths, linear_paths, faout, prefix):
    """
    Assemble cycle and linear paths with edge FASTA sequences and output to file.
    Function reconstructed from original make_fa_from_order, takes list parameters directly.
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
    n_seq = Seq("N" * 50)
    
    with open(faout, "w") as faout_stream:
        count = 0
        
        def write_paths(paths, tag):
            nonlocal count
            for path in paths:
                seq = ""
                for t in path:
                    if t == '':
                        continue
                    t = t.replace("ref", "")
                    node_name = t[0:-1]
                    orient = t[-1]
                    
                    if node_name not in record_dict:
                        print(f"Warning: Node '{node_name}' not found in {fain}", file=sys.stderr)
                        continue
                        
                    tmp_seq = record_dict[node_name].seq
                    if orient == '-':
                        tmp_seq = tmp_seq.reverse_complement()
                        
                    if seq == "":
                        seq = tmp_seq
                    else:
                        seq = seq + n_seq + tmp_seq
                        
                if seq != "":
                    count += 1
                    faout_stream.write(">" + prefix + "_phage_" + str(count) + "_" + tag + "\n")
                    faout_stream.write(str(seq) + "\n")

        # 保证 Cycle 排在前面，Linear 排在后面
        write_paths(cycle_paths, "cycle")
        write_paths(linear_paths, "linear")

def process_sample(path_file, graph_file, edge_fasta, out_fasta, prefix, trim_threshold, min_cycle_length):
    print(f"Checking path file: {path_file}", file=sys.stderr)
    print(f"Using graph file: {graph_file}", file=sys.stderr)
    print(f"Using edge fasta: {edge_fasta}", file=sys.stderr)
    print(f"Thresholds -> Max Trim: <= {trim_threshold} bp | Min Cycle Length: >= {min_cycle_length} bp", file=sys.stderr)

    if not os.path.exists(path_file):
        print(f"Error: Path file not found - {path_file}", file=sys.stderr)
        return
    if not os.path.exists(graph_file):
        print(f"Error: Graph file not found - {graph_file}", file=sys.stderr)
        return
    if not os.path.exists(edge_fasta):
        print(f"Error: Edge fasta file not found - {edge_fasta}", file=sys.stderr)
        return

    adjacency = parse_graph(graph_file)
    
    cycle_paths = []
    linear_paths = []

    with open(path_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "all" in line:
                continue

            path = [t for t in re.split(r'\s+', line) if t]
            is_circle, trimmed_path = is_circular_path_fuzzy(path, adjacency, trim_threshold, min_cycle_length)
            
            if is_circle:
                cycle_paths.append(trimmed_path)
            else:
                linear_paths.append(path)

    make_fa_from_paths(edge_fasta, cycle_paths, linear_paths, out_fasta, prefix)

    print(f"Total circular paths found: {len(cycle_paths)}", file=sys.stderr)
    print(f"Total linear paths found: {len(linear_paths)}", file=sys.stderr)
    print(f"Fasta successfully written to: {out_fasta}\n", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Check circular paths from assembly graph and generate fasta sequences."
    )
    
    parser.add_argument("path_file", help="Input path file containing contig orders")
    parser.add_argument("graph_file", help="Input assembly graph file (e.g., filtered_graph.txt)")
    parser.add_argument("edge_fasta", help="Input fasta file containing edge sequences")
    parser.add_argument("out_fasta", help="Output fasta file path")
    parser.add_argument("prefix", help="Prefix for the output sequence names (e.g., Sample_01)")
    
    parser.add_argument("--trim_threshold", type=int, default=300, 
                        help="Maximum length of contigs that can be trimmed from ends (default: 300)")
    parser.add_argument("--min_cycle_length", type=int, default=8000, 
                        help="Minimum physical length of the circular path (default: 8000)")

    args = parser.parse_args()

    process_sample(
        args.path_file, 
        args.graph_file, 
        args.edge_fasta, 
        args.out_fasta, 
        args.prefix, 
        args.trim_threshold, 
        args.min_cycle_length
    )

if __name__ == "__main__":
    main()