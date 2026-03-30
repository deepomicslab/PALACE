import argparse
import re
import subprocess
import tempfile
from collections import defaultdict

def find_index(lst, item):
    for idx, ele in enumerate(lst):
        if item == ele[2]:
            return ele[0]
    return -2


def create_sub_graph_from_remain_segs(removed_segs, full_segs):
    # remain_juncs = []
    # for junc in full_juncs:
    #     if junc not in removed_juncs:
            # remain_juncs.append(junc)
    tmp_remain_segs = []
    pure_segs = []
    for seg in full_segs:
        #print(full_segs[seg])
        is_add = True
        for rseg in removed_segs:
            if seg == rseg[1]:
                is_add = False
        if is_add:
            pure_segs.append([seg])
            tmp_remain_segs.append(f"SEG {seg} {' '.join(full_segs[seg])}")
    return pure_segs,tmp_remain_segs
def main(args):
    # Accessing the parsed arguments
    graph_file = args.graph
    prefix = args.prefix
    match_file = args.match
    samtools_command = args.samtools
    bam_file = args.bam
    blast_file = args.assembly_blast

    ref_pecent = parse_files(args.ref_percent_file)

    full_segs, full_juncs = parse_graph_file(graph_file)
    #seg_depths = get_segs_depth(full_segs.keys(), bam_file, samtools_command)
    seg_depths = pysam.TabixFile(bam_file)
    graph_dict, similar_refs = parse_match_file(match_file, ref_pecent)
    ref_order = parse_blast(blast_file)
    
    with open(args.ref_same_out, "w") as f:
        # Sort dictionary keys to ensure consistent output order
        for key in sorted(similar_refs.keys()):
            item = similar_refs[key]
            f.write(",".join(item))
            f.write("\n")
            
    # Extract similar_refs values in sorted order
    similar_refs_list = [item for key in sorted(similar_refs.keys()) for item in similar_refs[key]]
    added_segs = []
    orders = []
    
    # Traverse graph_dict in sorted key order
    for ref_key, ref_segs in sorted(graph_dict.items()):
        if ref_key not in similar_refs_list:
            continue
        if ref_key in ref_order.keys():
            orders = ref_order[ref_key]
        updated_ref_segs = update_segs_with_depth(ref_segs, seg_depths, full_segs)
        if len(updated_ref_segs) == 0:
            continue
            
        with open(prefix +"_ref"+ ref_key+"ref.second", "w") as f:
            ref_juncs = parse_juncs_from_segs(ref_segs, full_juncs)
            for seg in updated_ref_segs:
                added_segs.append(seg)
                seg_order = find_index(orders, seg[1])
                if seg_order == -2:
                    seg[-1] = '-1'
                f.write(" ".join(seg)+" "+str(seg_order)+"\n")
            
            # Convert set to sorted list to eliminate hash randomization
            for junc in sorted(ref_juncs):
                f.write(junc+"\n")
                
    pure_segs, remain_segs = create_sub_graph_from_remain_segs(added_segs, full_segs)
    
    with open(prefix +"_refremain"+"ref.second", "w") as f:
        ref_juncs = parse_juncs_from_segs(pure_segs, full_juncs)
        for seg in remain_segs:
            added_segs.append(seg)
            f.write(seg+" -1"+"\n")
            
        # Convert set to sorted list for output
        for junc in sorted(ref_juncs):
            f.write(junc+"\n")


def get_segs_depth(segs, bam, samtools):
    """
    Calculate depths for multiple segments using a BED file for efficiency.
    """
    seg_depths = {}

    # Collect all contigs to query
    contigs = [item for item in segs]

    # Get depths for all contigs via a single samtools depth call
    depth_data = run_samtools_depth_with_bed(bam, contigs)
    print("depth donw")

    if depth_data:
        for contig in contigs:
            if contig in depth_data:
                depths = depth_data[contig]
                average_depth = calculate_average_depth(depths)
                seg_depths[contig] = (average_depth, len(depths))
            else:
                # Handle case where no depth data is returned for a contig
                spade_contig = contig.split("_")
                seg_depths[contig] = (spade_contig[-1], spade_contig[-3])

    return seg_depths

def parse_juncs_from_segs(segs, full_juncs):
    ref_juncs = set()
    segs_1d = [item for row in segs for item in row]
    for jun_key, jun_str in full_juncs.items():
        if jun_key[0] in segs_1d and jun_key[2] in segs_1d:
            ref_juncs.add(" ".join(jun_str))
    return ref_juncs


import pysam

def run_samtools_depth_with_bed(depth_gz_file, contigs):
    """
    Parse a tabix-indexed samtools depth.gz file for specific contigs.

    This function reads a tabix-indexed depth.gz file and extracts depth information for the specified contigs.

    Args:
        depth_gz_file (str): Path to the tabix-indexed depth.gz file (compressed output from samtools depth).
        contigs (list): List of contig names (e.g., ["chr1", "chr2"]).

    Returns:
        dict: A dictionary where keys are contigs, and values are lists of depth values for the specified contigs.
    """
    # Dictionary to store the extracted depth values
    depth_dict = {}

    # Open the tabix-indexed depth.gz file using pysam
    try:
        with pysam.TabixFile(depth_gz_file) as tabix:
            for contig in contigs:
                # Check if the contig exists in the depth.gz file
                if contig in tabix.contigs:
                    # Initialize an empty list for the contig
                    depth_dict[contig] = []

                    # Fetch all records for the contig
                    for record in tabix.fetch(contig):
                        # Parse the record (format: contig, position, depth)
                        _, _, depth = record.split("\t")
                        depth_dict[contig].append(int(depth))
                else:
                    print(f"Contig {contig} not found in the depth.gz file.")
    except Exception as e:
        print(f"Error reading the depth.gz file: {e}")

    return depth_dict




def calculate_average_depth(depths):
    total_depth = sum(depths)
    if len(depths) == 0:
        return 0
    average_depth = total_depth / len(depths)
    return average_depth

import pysam

def update_segs_with_depth(segs, depth_tabix, seg_gene_scores):
    """
    Update segments with depth information using a Tabix-indexed depth.gz file.
    If no depth information is found for a contig, fallback to extracting depth and length from contig name.

    Args:
        segs (list): List of contig names (e.g., ["contig_1_500_30", "contig_2_600_20"]).
        depth_tabix (pysam.TabixFile): Tabix-indexed depth.gz file for querying depth information.
        seg_gene_scores (dict): Dictionary with contig as key and gene scores as value.

    Returns:
        list: Updated segments with copy number and other information.
    """
    total_depths = 0
    total_lens = 0
    seg_depths = {}  # Store calculated depths for each contig
    final_segs = []

    # Step 1: Calculate depths for all contigs and accumulate total depths and lengths
    for item in segs:
        contig = item[0]
        depths = []
        contig_len = 0
        avg_depth = 0

        # Try fetching depth information for the contig
        try:
            # Fetch all depth values for the contig from the tabix file
            for record in depth_tabix.fetch(contig):
                _, _, depth = record.split("\t")
                depths.append(int(depth))
        except ValueError:
            # If no depth information is found, fall back to parsing the contig name
            spade_contig = contig.split("_")
            avg_depth = float(spade_contig[-1])  # Extract average depth from contig name
            contig_len = int(spade_contig[-3])  # Extract contig length from contig name
            seg_depths[contig] = (avg_depth, contig_len)
            total_depths += avg_depth * contig_len
            total_lens += contig_len
            continue

        # Calculate the average depth and total base count for the contig
        if depths:
            avg_depth = sum(depths) / len(depths)
            contig_len = len(depths)
            seg_depths[contig] = (avg_depth, contig_len)
            total_depths += avg_depth * contig_len
            total_lens += contig_len

    # If total length is 0, return an empty list
    if total_lens == 0:
        return []

    # Step 2: Calculate the total average depth across all contigs
    total_average_depth = total_depths / total_lens

    # Step 3: Use the precomputed depths for each contig to compute copy number and update segments
    for item in segs:
        contig = item[0]
        # Use precomputed contig depth information
        if contig in seg_depths:
            avg_depth, contig_len = seg_depths[contig]
            copy_num = round(avg_depth / total_average_depth)
            if copy_num == 0:
                copy_num = 1  # Ensure copy number is at least 1

            # Append the updated segment information
            final_segs.append([
                "SEG",
                contig,
                str(avg_depth),
                str(copy_num),
                seg_gene_scores[contig][2] if contig in seg_gene_scores else "0",
                seg_gene_scores[contig][3] if contig in seg_gene_scores else "0",
                "1"
            ])

    return final_segs


def parse_graph_file(file_path):
    segs = {}
    juncs = {}

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if not parts:
                continue  # skip empty lines

            if parts[0] == "SEG":
                # Assuming we want to store the whole line for SEGs,
                # you could also extract specific fields if needed
                segs[parts[1]] = parts[2:]
            elif parts[0] == "JUNC":
                juncs[(parts[1], parts[2], parts[3], parts[4])] = parts
                # Assuming we want to store the whole line for JUNCs,
                # you could also extract specific fields if needed
                # juncs.append(parts[1:])
    return segs, juncs
def parse_match_file(file_path, ref_pecent):
    similar_refs = {}
    graph_dict = {}

    # Regular expression to find edges with their orientations
    edge_pattern = re.compile(r'(EDGE_[\w_]+_cov_[\d.]+)([+-])')

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if not parts:
                continue  # Skip empty lines
            # The last part of the line is the sequence identifier,
            # everything before that are concatenated edges with orientation
            seq_id = parts[-1]
            #if float(ref_pecent[seq_id]) < 0.9:
            #    continue
            if parts[0] in similar_refs.keys():
                similar_refs[parts[0]].append(parts[-1])
            else:
                similar_refs[parts[0]] = [parts[-1]]
            edge_string = ' '.join(parts[:-1])  # Join all but last part into a single string

            # Extract edges with orientations using regex
            edges_with_orientation = [(match.group(1), match.group(2)) for match in edge_pattern.finditer(edge_string)]

            # Add to dictionary, initializing if not already present
            if seq_id not in graph_dict:
                graph_dict[seq_id] = []
            graph_dict[seq_id].extend(edges_with_orientation)
    for k,refs in similar_refs.items():
        max_percent = 0;
        max_percent_ref = ''
        for ref in refs[:]:
            if max_percent < ref_pecent[ref]:
                max_percent = ref_pecent[ref]
                max_percent_ref = ref
            else:
                if ref_pecent[ref] < 0.85:
                    similar_refs[k].remove(ref)
        if len(similar_refs[k]) == 0:
            similar_refs[k].append(max_percent_ref)

    return graph_dict, similar_refs

def parse_blast(blast_file):
    # Dictionary to hold the reference and their corresponding query contigs
    reference_dict = defaultdict(list)

    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue  # Skip incomplete lines

            query_id = parts[0]
            subject_id = parts[1]
            s_start = min(int(parts[8]), int(parts[9]))
            s_end = max(int(parts[8]), int(parts[9]))
            sublen = int(parts[13])
            querylen = int(parts[12])
            found = False
            current_len = s_end - s_start
            for idx, item in enumerate(reference_dict[subject_id]):
                if query_id == item[2]:
                    if abs(s_start - s_end) > abs(item[0] - item[1]):
                        reference_dict[subject_id][idx] = (s_start, s_end, query_id, reference_dict[subject_id][idx][3] + current_len/querylen)
                    elif s_start - 1 < 10:
                        if sublen - item[1] < 50: # circular
                            if s_end == int(parts[9]):
                                reference_dict[subject_id][idx] = (0, s_end, query_id, reference_dict[subject_id][idx][3] + current_len/querylen)
                            else:
                                reference_dict[subject_id][idx] = (-1, s_end, query_id, reference_dict[subject_id][idx][3] + current_len/querylen)
                    else:
                        reference_dict[subject_id][idx] = (reference_dict[subject_id][idx][0], reference_dict[subject_id][idx][1], reference_dict[subject_id][idx][2], reference_dict[subject_id][idx][3] + current_len/querylen)
                    found = True
            if not found:
            # Store the query id with its start position in the reference
                reference_dict[subject_id].append((s_start, s_end, query_id, current_len/querylen))

    # Sort the query contigs for each reference based on the start position in the reference
    updated_data = {
    key: [(-2, b, c, d) if d < 0.5 else (a, b, c, d) for (a, b, c, d) in value]
    for key, value in reference_dict.items()
    }
    for subject_id in updated_data:
        updated_data[subject_id].sort()
    # Output the sorted order of query contigs for each reference
    # for subject_id, alignments in reference_dict.items():
    #     print(f"Reference: {subject_id}")
    #     for alignment in alignments:
    #         s_start, s_end, query_id = alignment
    #         print(f"Query contig: {query_id}, Position: {s_start}-{s_end}")
    return updated_data
def parse_files(data_file_path):
    contig_dict = {}

    with open(data_file_path, 'r') as data_file:
        data_lines = data_file.readlines()

        for data_line in data_lines:
            line_array = data_line.split('\t')
            contig_name = line_array[0]
            float_value = float(line_array[-1])
            contig_dict[contig_name] = float_value

    return contig_dict
if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process some files and a command.")

    # Adding arguments
    parser.add_argument("graph", help="The file path to the graph file.")
    parser.add_argument("prefix", help="The prefix for processing or output.")
    parser.add_argument("match", help="The file path to the match file.")
    parser.add_argument("samtools", help="The samtools command or path.")
    parser.add_argument("bam", help="bam file path.")
    parser.add_argument("assembly_blast", help="assembly blast file path.")
    parser.add_argument("ref_same_out", help="similar refs.")
    parser.add_argument("ref_percent_file", help="refs percent")

    # Parse arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
