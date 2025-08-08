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
        for item in similar_refs.values():
            f.write(",".join(item))
            f.write("\n")
    similar_refs_list = [item for sublist in similar_refs.values() for item in sublist]
    added_segs = []
    orders = []
    for ref_key, ref_segs in graph_dict.items():
        if ref_key not in similar_refs_list:
            continue
        if ref_key in ref_order.keys():
            orders = ref_order[ref_key]
        updated_ref_segs = update_segs_with_depth(ref_segs, seg_depths, full_segs)
        if len(updated_ref_segs) == 0:
            continue
        with open(prefix +"_ref"+ ref_key+"ref.second", "w") as f:
            # print(ref_key, ref_segs,"------", updated_ref_segs)
            # print("================")
            ref_juncs = parse_juncs_from_segs(ref_segs, full_juncs)
            for seg in updated_ref_segs:
                added_segs.append(seg)
                seg_order = find_index(orders, seg[1])
                if seg_order == -2:
                    seg[-1] = '-1'
                f.write(" ".join(seg)+" "+str(seg_order)+"\n")
            for junc in ref_juncs:
                # print(junc)
                f.write(junc+"\n")
    pure_segs,remain_segs = create_sub_graph_from_remain_segs(added_segs, full_segs)
    with open(prefix +"_refremain"+"ref.second", "w") as f:
        #print(remain_segs)
        #updated_ref_segs = update_segs_with_depth(remain_segs, seg_depths, full_segs)
        ref_juncs = parse_juncs_from_segs(pure_segs, full_juncs)
        for seg in remain_segs:
            added_segs.append(seg)
            #seg_order = find_index(orders, seg[1])
            #if seg_order == -2:
            #    seg[-1] = '-1'
            f.write(seg+" -1"+"\n")
        for junc in ref_juncs:
            # print(junc)
            f.write(junc+"\n")
#def get_segs_depth(segs, bam, samtools):
#    seg_depths = {}
#    total_depths = []
#    for item in segs:
#        contig = item
#        depths = run_samtools_depth(bam, contig, samtools)
#        if depths:
#            average_depth = calculate_average_depth(depths)
#            seg_depths[contig] = (average_depth, len(depths))
#        else:
#            spade_contig = contig.split("_")
#            seg_depths[contig] = (spade_contig[-1], spade_contig[-3])
#    return seg_depths

    # depths = run_samtools_depth(bam, segs, samtools)
    # seg_depths = dict(zip(segs, depths))
    # return seg_depths
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

# def parse_juncs_from_segs(segs, full_juncs):
#     ref_juncs = set()
#     for item in segs:
#         for jun_key, jun_str in full_juncs.items():
#             if item[0] in jun_key:
#                 ref_juncs.add(" ".join(jun_str))
#     return ref_juncs
#def run_samtools_depth(input_bam, region, samtools):
#    cmd = f'{samtools} depth -r {region} {input_bam}'
#    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#
#    if result.returncode != 0:
#        print(f"Error running samtools depth: {result.stderr}")
#        return None
#
#    depths = [int(line.split('\t')[2]) for line in result.stdout.splitlines()]
#    return depths

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


#import gzip

#def run_samtools_depth_with_bed(depth_gz_file, contigs):
#    """
#    Parse a samtools depth.gz file for specific contigs.
#
#    This function reads an indexed depth.gz file and extracts depth information for the specified contigs.
#
#    Args:
#        depth_gz_file (str): Path to the depth.gz file (compressed output from samtools depth).
#        contigs (list): List of contig names (e.g., ["chr1", "chr2"]).
#
#    Returns:
#        dict: A dictionary where keys are contigs, and values are lists of depth values for the specified contigs.
#    """
#    # Convert contigs list to a set for faster lookups
#    contig_set = set(contigs)
#
#    # Dictionary to store the extracted depth values
#    depth_dict = {}
#
#    # Open and read the depth.gz file
#    with gzip.open(depth_gz_file, "rt") as f:
#        for line in f:
#            contig, pos, depth = line.split("\t")  # Each line contains contig, position, depth
#            pos, depth = int(pos), int(depth)
#
#            # Check if the contig is in the specified list of contigs
#            if contig in contig_set:
#                if contig not in depth_dict:
#                    depth_dict[contig] = []
#                depth_dict[contig].append(depth)
#
#    return depth_dict

#def run_samtools_depth_with_bed(input_bam, regions, samtools):
#    """
#    Run samtools depth using a BED file for specifying regions.
#    This function creates a BED file named based on the input BAM file, runs samtools depth, and removes the BED file afterward.
#
#    Args:
#        input_bam (str): Path to the input BAM file.
#        regions (list): List of regions in BED format (e.g., ["contig1\t0\t1000", "contig2\t500\t1500"]).
#        samtools (str): Path to the samtools executable.
#
#    Returns:
#        dict: A dictionary where keys are contigs, and values are lists of depth values for the specified regions.
#    """
#    # Create a BED file named based on the input BAM file
#    bed_file_path = f"{input_bam}.bed"
#    try:
#        # Write regions to the BED file
#        with open(bed_file_path, "w") as bed_file:
#            for region in regions:
#                bed_file.write(f"{region}\n")  # Write regions in BED format
#            bed_file.flush()  # Flush the buffer to ensure data is written to disk
#            os.fsync(bed_file.fileno())  # Ensure the file is physically written to disk
#
#        # Run samtools depth with the BED file
#        cmd = f'{samtools} depth -b {bed_file_path} {input_bam}'
#        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#
#        if result.returncode != 0:
#            print(f"Error running samtools depth: {result.stderr}")
#            return None
#
#        # Parse the output into a dictionary
#        depth_dict = {}
#        for line in result.stdout.splitlines():
#            contig, pos, depth = line.split('\t')  # Extract contig, position, and depth
#            depth = int(depth)
#            if contig not in depth_dict:
#                depth_dict[contig] = []
#            depth_dict[contig].append(depth)
#
#        return depth_dict
#    finally:
#        pass
        # Clean up the BED file
        #if os.path.exists(bed_file_path):
        #    os.remove(bed_file_path)

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
#def update_segs_with_depth(segs, seg_depths, seg_gene_scores):
#    total_depths = 0
#    total_lens = 0
#    average_contig_depth = {}
#    final_segs = []
#    for item in segs:
#        contig = item[0]
#        avg_depth, contig_len = seg_depths[contig]
#        total_depths = total_depths + float(avg_depth) * int(contig_len)
#        total_lens = total_lens + int(contig_len)
#    if total_lens == 0:
#        return []
#    total_average_depth = total_depths/total_lens
#    # print(segs, total_average_depth)
#    for k in seg_depths.keys():
#        in_segs = False
#        for seg in segs:
#            if k in seg:
#                in_segs = True
#        if not in_segs:
#            continue
#        i_arr = k.split()
#        i_depth = seg_depths[k]
#        copy_num = round(float(i_depth[0])/total_average_depth)
#        if copy_num == 0:
#            copy_num = 1
#        # final_segs.append("SEG " + k + " " +k.split('_')[-1] + " "+str(copy_num) + " 0 0 1")
#        final_segs.append(["SEG",k, str(i_depth[0]), str(copy_num) ,seg_gene_scores[k][0], seg_gene_scores[k][1], "1"])
#    return final_segs

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
