import pysam
import sys
import re

def round_with_threshold(num, threshold=0.5):
    if num - int(num) < threshold:
        if int(num) == 0:
            return 1
        return int(num)
    else:
        return int(num) + 1

def average_depth(samfile, contig_name):
    total_depth = 0
    length = 0

    for pileupcolumn in samfile.pileup(contig_name):
        total_depth += pileupcolumn.n
        length += 1
    # Calculate and return average depth
    if length > 0:
        return total_depth / length
    else:
        return 0

def parse_pair(p1,p2):
    p1_strand = p1[2]
    p2_strand = p2[2]
    result = None
    if p1_strand == p2_strand:
        if (p1[0] == 0 and p2[0] == 0) or (p1[0] > 0 and p2[0] > 0):
            return result
        if p1[0] == 0 and p2[0] > 0:
            if p1[3] < 100 and p2[4] < 100:
                result = [p2[1], "+", p1[1], "+"]
        elif p1[0] > 0 and p2[0] == 0:
            if p1[4] < 100 and p2[3] < 100:
                result = [p1[1], "+", p2[1], "+"]
    else:
        if (p1[0] == 0 and p2[0] > 0) or (p1[0] > 0 and p2[0] == 0):
            return result
        # if p1[0] == 0 and p2[0] == 0:
        #     if p1[3]
        result = [p1[1], "+", p2[1], "-"]
    if result is not None:
        return " ".join(result)
    else:
        return None
def generate_one_side(query_name, query_start, query_end, ref_name, ref_length, cigar_string, strand):
    max_soft_index = parse_cigar(cigar_string)
    if max_soft_index == -1:
        return None
    # max_soft_index = parse_cigar(cigar_string)
    # print(strand, max_soft_index, query_start, query_end)
    return (max_soft_index, ref_name, "+", query_start, ref_length - query_end)
    # if strand == "+":
    #     if max_soft_index == 0:
    #         if query_start < 100:
    #             return (0, ref_name, "+", query_start, query_end)
    #         else:
    #             return (None, None, None)
    #     elif max_soft_index > 0:
    #         if ref_length - query_end < 100:
    #             return (1, ref_name, "+")
    #         else:
    #             return (None, None, None)
    # else:
    #     if max_soft_index == 0:
    #         if query_start < 100:
    #             return (1, ref_name, "-")
    #         else:
    #             return (None, None, None)
    #     elif max_soft_index > 0:
    #         if ref_length - query_end < 100:
    #             return (0, ref_name, "+")
    #         else:
    #             return (None, None, None)
    # return (None, None, None)

def parse_cigar(cigar_string):
    # Regular expression to match a number followed by a letter
    cigar_tuples = re.findall("(\d+)([MIDNSHP=X])", cigar_string)
    cigar_values = [int(value) for value, key in cigar_tuples]
    cigar_keys = [key for value, key in cigar_tuples]

    # Check for soft clips and update max_soft_clip_len if a longer soft clip is found
    max_soft_clip_len = -1
    for i in range(len(cigar_keys)):
        if cigar_keys[i] == "S":
            max_soft_clip_len = max(max_soft_clip_len, cigar_values[i])
    if max_soft_clip_len != -1:
        max_soft_index = cigar_values.index(max_soft_clip_len)
        return max_soft_index
    else:
        return -1

# Open the SAM file
samfile = pysam.AlignmentFile(sys.argv[1],"r")
total_avg_depth = int(float(sys.argv[2]))
# Iterate over each record in the SAM file
read_side = {}
ref_lengths = samfile.lengths
for read in samfile:
    # Skip if not a supplementary alignment
    if read.is_secondary:
        continue
    # Get the query name, reference name, CIGAR string, and strand
    query_name = read.query_name
    query_start = read.reference_start
    query_end = read.reference_end
    ref_length = ref_lengths[read.reference_id]
    ref_name = samfile.get_reference_name(read.reference_id)
    #cigar_string = read.cigarstring
    cigar_string = read.cigarstring.replace("H", "S")
    strand = '-' if read.is_reverse else '+'
    if not read.has_tag('SA'):
        continue
    sa_tag = read.get_tag("SA")
    if sa_tag is None or len(sa_tag.split(";")) > 2:
        continue
    side_info = generate_one_side(query_name, query_start, query_end, ref_name, ref_length, cigar_string, strand)
    if query_name not in read_side.keys():
        read_side[query_name] = []
    read_side[query_name].append(side_info)
    # sa_array = sa_tag.split(":")[2].split(",")
    # supp_ref_name = sa_array[0]
    # supp_query_start = int(sa_array[1])
    # supp_strand = sa_array[2]
    # supp_cigar = sa_array[3]
    # supp_cigar_dict = parse_cigar(supp_cigar)


    # Print the information
    # print(f"Query Name: {query_name}, Reference Name: {ref_name}, CIGAR String: {cigar_string}, Strand: {strand}")

juncs = {}
for key, value in read_side.items():
    if value[0] == None or value[1] == None:
        continue
    junc = parse_pair(value[0], value[1])
    if junc == None:
        continue
    if junc in juncs.keys():
        juncs[junc] += 1
    else:
        juncs[junc] = 1
contig_names = samfile.references
for key in contig_names:
    avg_depth = average_depth(samfile, key)
    copy1 = float(avg_depth)/float(total_avg_depth)
    copy_num = round_with_threshold(copy1, 0.7)
    print('SEG '+key+' '+str(avg_depth)+' '+str(copy_num))
for key,value in juncs.items():
    print('JUNC '+key+" "+ str(value))

# Close the SAM file
samfile.close()
