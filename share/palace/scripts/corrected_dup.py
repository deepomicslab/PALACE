import sys
from pyfaidx import Fasta
import re
from Bio import SeqIO
import re
import sys
from Bio.Seq import Seq
import re
import subprocess
import os
import copy
from collections import Counter
from itertools import chain
def get_path_len(path):
    sum = 0
    for item in path:
        if item.startswith("EDGE"):
            arr = item.split("_")
            sum +=int(arr[3])
    return sum
def split_list_on_element(input_list, A):
    # Find the indices of element A in the list
    indices = [i for i, elem in enumerate(input_list) if A in elem]

    # Add the last index for easy slicing
    indices.append(len(input_list))

    # Create sublists by slicing the list using the indices
    sublists = [input_list[indices[i]:indices[i + 1]] for i in range(len(indices) - 1)]

    # Convert sublists to tuples and count their occurrences
    sublist_counts = Counter(tuple(sublist) for sublist in sublists)
    return sublist_counts

def concatenate_sublists_with_counts(sublist_counts):
    repeated_sublists = [list(sublist) * count for sublist, count in sublist_counts.items()]
    flattened_list = list(chain.from_iterable(repeated_sublists))

    return flattened_list

def merge_repeat(lst):
    # Count the occurrences of characters in the string, if we need?
    s = [item.replace("-","").replace("+","") for item in lst]
    # s = lst
    item_counts = Counter(s)

    # Find the most frequent item
    most_frequent_item = max(item_counts, key=item_counts.get)

    # Find the first occurrence of the most frequent item
    item_index = s.index(most_frequent_item)

    # Reformat the string to start with the most frequent item: cabaea to acabae
    reformatted_list = lst[item_index:] + lst[:item_index]
    sublist_counts = split_list_on_element(reformatted_list, most_frequent_item)
    return concatenate_sublists_with_counts(sublist_counts)



def run_samtools_depth(input_bam, region):
    cmd = f'samtools depth -r {region} {input_bam}'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        #print(f"Error running samtools depth: {result.stderr}")
        return None

    depths = [int(line.split('\t')[2]) for line in result.stdout.splitlines()]
    return depths

def calculate_average_depth(depths):
    total_depth = sum(depths)
    if len(depths) == 0:
        return 0
    average_depth = total_depth / len(depths)
    return average_depth

def get_min_copy_seg(unit_seg, seg_copies):
    # segs = unit_seg.split()
    min_seg = ""
    min_copy = 10000
    for item in unit_seg:
        item = item.replace("+","").replace("-","")
        if item not in seg_copies.keys():
            copy = 1
        else:
            copy = seg_copies[item]
        if copy < min_copy:
            min_seg = item
            min_copy = copy
    return min_seg, min_copy

def non_dup_item(ori_arr, unit_cycles):
    ori_arr_str = "\t".join(ori_arr).replace("+","").replace("-","")
    unit_cycle_str = ["\t".join(item).replace("+","").replace("-","") for item in unit_cycles]
    for item in unit_cycle_str:
        ori_arr_str.replace(item,"")
    return ori_arr_str.split("\t")

def calculate_real_copy_for_cycle(unit_seg, seg_copies, non_unit_part):
    min_seg, min_copy = get_min_copy_seg(unit_seg, seg_copies)
    other_count = non_unit_part.count(min_seg)
    real_copy = min_copy - other_count
    if real_copy < 1:
        real_copy = 1
    return real_copy

def get_depth(all_non_dup_segs,unit_cycle_segs, non_unit_part, bam, first_item):
    unit_copies = []
    seg_len_depth = {}
    total_depths = []
    total_positions = 0
    average_contig_depth = {}
    final_segs = []
    seg_depth_dict = {}
    for item in all_non_dup_segs:
        contig = item.replace('-','').replace('+','')
        depths = run_samtools_depth(bam, contig)
        if depths:
            average_depth = calculate_average_depth(depths)
            seg_len_depth[contig] = [average_depth, len(depths)]
            total_depths.extend(depths)
            total_positions += len(depths)
    total_average_depth = calculate_average_depth(total_depths)
    for k in seg_len_depth.keys():
        seg_depth_dict[k] = round(seg_len_depth[k][0]/total_average_depth)
    for unit_seg in unit_cycle_segs:
        seg_unit_copy = calculate_real_copy_for_cycle(unit_seg, seg_depth_dict, non_unit_part)
        unit_copy = round(seg_unit_copy)
        if unit_copy == 0:
            unit_copy = 1
        unit_copies.append(unit_copy)
    k = first_item.replace('-','').replace('+','')
    return_depth=0
    if k in seg_depth_dict.keys():
        return_depth = seg_depth_dict[first_item.replace('-','').replace('+','')]
    return unit_copies, return_depth


def reformat_cycle(s):
    ori_s = copy.deepcopy(s)
    n = len(s)
    longest_substring_index = -1

    for i in range(n // 2 + 1):
        if s[:i] == s[-i:]:
            longest_substring_index = i
    #print(longest_substring_index)
    if longest_substring_index != -1:
        return s[len(s) - longest_substring_index:] + s[0: len(s) - longest_substring_index]
    if ori_s == s:
        s = merge_repeat(ori_s)
    return s

def are_cyclically_equal(s1, s2):
    if s1 in s2:
        return True
    #if len(s1) != len(s2):
    #    return False
    concatenated = s1 +"\t"+ s1
    return s2 in concatenated

def find_consecutive_repeats(s, min_repeat_len=2):
    repeats = set()
    for repeat_len in range(1, len(s) // 2 + 1):
        for start in range(0, len(s) - repeat_len * 2 + 1):
            found = False
            count = 1
            while s[start:start + repeat_len] == s[start + repeat_len * count:start + repeat_len * (count + 1)]:
                found = True
                count += 1
            tag = True
            if found and count >= min_repeat_len:
                for item in repeats:
                    if are_cyclically_equal(item,"\t".join(s[start:start + repeat_len])):
                        tag = False
                        break
                if tag:
                    repeats.add("\t".join(s[start:start + repeat_len]))
    return [item.split("\t") for item in repeats]

def find_consecutive_repeats_rotated(s, min_repeat_len=2):
    n = len(a)
    final_list = a
    repeats = set()
    for i in range(n):
        rotated_lst = a[i:] + a[:i]
        current_repeats =remove_dup_record(rotated_lst)
        if len(current_list) < len(final_list):
            final_list = current_list
    return final_list
#split result into fasta

# arg 1: fasta
# arg 2: result.txt

# split fasta with result.txt
# def choose_faidx_key(lst, edge_fasta):
#     lst_len = sum([int(item.split("_")[3]) for item in lst])
#     for edge_key in edge_fasta.keys():
#         if len(edge_fasta[edge_key]) == lst_len:
#             return edge_key
#     return None
def get_prefix(edge_fasta):
    return list(edge_fasta.keys())[0].split("_")[0]

def split_fasta(edge_fasta_file, final_all_file, edge_fasta_out_file):
    edge_fasta = Fasta(edge_fasta_file)
    edge_txt = open(final_all_file)
    edge_out = open(edge_fasta_out_file, "w")
    writed_edge = []
    idx = 1
    prefix = get_prefix(edge_fasta) # SAMEA728574_phage_3, prefix:SAMEA728574
    for line in edge_txt.readlines():
        prev_len = 0
        line_arr = re.split(r'\s+', line.strip())
        # # #print(line_arr)
        edge_fasta_key = prefix+"_phage_" + str(idx)
        # #print(edge_fasta_key, len(edge_fasta[edge_fasta_key]))
        for item in line_arr:
            item_dir = item[-1]
            item_len = get_seg_len(item)
            item_fa = edge_fasta[edge_fasta_key][prev_len:item_len + prev_len]
            if item_dir == "-":
                item_fa_rev = item_fa.reverse.complement
            prev_len = item_len + prev_len
            if item[:-1] in writed_edge:
                continue
            edge_out.write(">"+item[:-1]+"\n")
            writed_edge.append(item[:-1])
            if item_dir == "-":
                edge_out.write(str(item_fa_rev) + "\n")
            else:
                edge_out.write(str(item_fa) + "\n")
        idx = idx + 1
    edge_txt.close()
    edge_out.close()
    return prefix

# def remove_dup(edge_txt_file):
#     line_lst = []
#     for line in open(edge_txt_file).readlines():
#         sub_lst = []
#         line_arr = re.split(r"\s+", line)
#         for item in line_arr:
#             sub_lst.append(item[:-1])

def find_sublist_indexes(A, B):
    if not A or not B:
        return -1, -1

    first_index = -1
    last_index = -1

    for i in range(len(B) - len(A) + 1):
        if B[i:i + len(A)] == A:
            if first_index == -1:
                first_index = i
            last_index = i

    return first_index, last_index + len(A)

def count_item_ignor_direction(lst,ele):
    count = 0
    ele = ele.replace("+","").replace("-","")
    for item in lst:
        if ele in item:
            count=count+1
    return count
def get_contig_len_for_arr(lst):
    final_len = 0
    for item in lst:
        final_len = final_len + get_seg_len(item)
    return final_len
def push_back_cycle_copies(unit_cycles, unit_copies, line_arr, first_item_copy):
    for i in range(len(unit_cycles)):
        unit_item = unit_cycles[i] + unit_cycles[i]
        unit_copy = unit_copies[i]
        # at least one copy
        if unit_copy < 1:
            unit_copy = 1
        start_idx, end_idx = find_sublist_indexes(unit_item, line_arr)
        line_arr = line_arr[:start_idx] + unit_cycles[i] * unit_copy + line_arr[end_idx:]
    # check copy split cycle
    first_item_count_in_liner_arr = count_item_ignor_direction(line_arr, line_arr[0])
    if abs(first_item_count_in_liner_arr - first_item_copy) <= 1:
        return line_arr
    # if line_arr contains larger copy
    sublist_counts = split_list_on_element(line_arr, line_arr[0])
    final_list = []
    final_list_len = 0
    for sublist,count in sublist_counts.items():
        current_len = get_contig_len_for_arr(sublist)
        if current_len > final_list_len:
            final_list = sublist
            final_list_len = current_len
    return final_list

# filter cycle result
def filter_cycle(cycle_file, cycle_out_file,bam):
    cycle_out = open(cycle_out_file, "w")
    tmp_result = []
    line_count = 0
    ori_cycle_result = []
    for line in open(cycle_file).readlines():
        line_count = line_count + 1
        # sub_lst = []
        line_arr = re.split(r"\s+", line.strip())
        ori_cycle_result.append(line_arr)
        # for each record, filter cycle result like : ABABABAC to ABAC
#1. find unit cycle first, and calculate the non dup element segment depth and then, calculate the unit cycle copy, and correct the copy
        line_arr = reformat_cycle(line_arr)
        # the first item is the most occurrences  element
        first_item = line_arr[0]
        unit_cycles = find_consecutive_repeats(line_arr)
        #print("cycles",unit_cycles)
        non_unit_part = non_dup_item(line_arr, unit_cycles)
        unit_copies, first_item_copy = get_depth(set(line_arr),unit_cycles, non_unit_part ,bam, first_item)
        corrected_line_arr = push_back_cycle_copies(unit_cycles, unit_copies, line_arr, first_item_copy)
        tmp_result.append(corrected_line_arr)
    keeped_idx = set(range(0, len(tmp_result)))
    for i in range(0, len(tmp_result)):
        if i not in keeped_idx:
            continue
        for j in range(i, len(tmp_result)):
            if i == j:
                continue
            if j not in keeped_idx:
                continue
            similar, idx = is_similar(tmp_result[i], tmp_result[j])
            if similar:
                if idx == 0:
                    keeped_idx.remove(j)
                else:
                    keeped_idx.remove(i)
                    break
    final_result = [tmp_result[i] for i in keeped_idx]
    for row in final_result:
        line = ' '.join(str(elem) for elem in row) + '\n'
        cycle_out.write(line)
    cycle_out.close()
    return line_count, final_result, ori_cycle_result

# if lst1 contains lst2 ignore the order
def is_same(lst1, lst2):
    return set(lst2).issubset(set(lst1))

def is_similar(lst1, lst2):
    # # #print(lst1, lst2, "lst1 lst2")
    # lst1_sorted = sorted(lst1, reverse=True)
    lst1_lens = [get_seg_len(item) for item in lst1]

    # lst2_sorted = sorted(lst2, reverse=True)
    lst2_lens = [get_seg_len(item) for item in lst2]

    lst1_lens_sum = sum(set(lst1_lens))
    lst2_lens_sum = sum(set(lst2_lens))
    # lst1_lens_sum = sum(lst1_lens)
    # lst2_lens_sum = sum(lst2_lens)
    intersections = set(lst1_lens).intersection(lst2_lens)

    if sum(intersections)/lst1_lens_sum >= 0.9 or sum(intersections)/lst2_lens_sum >= 0.9:
        # if too long, choose the shorter one. if short, choose the longer one
        # if min(lst1_lens_sum, lst2_lens_sum) > 100000:
        #     return
        if lst1_lens_sum > lst2_lens_sum:
            return True, 0
        else:
            return True, 1
    return False, -1


def remove_dup_record2(lst):
    input_list = lst + [lst[0]]
    cycles = []
    final_cycles = []
    n = len(input_list)
    first_element = input_list[0]

    for i in range(n):
        if input_list[i] == first_element:
            for j in range(i + 2, n + 1):
                if input_list[j - 1] == first_element:
                    cycle = input_list[i:j]
                    cycles.append(cycle)
    cycles = sorted(cycles, key=len, reverse=True)
    unit_cycle_idx = set(range(0, len(cycles)))
    for i in range(0,len(cycles)):
        c1 = cycles[i]
        for j in range(i,len(cycles)):
            if j == i:
                continue
            c2 = cycles[j]
            if contains_sublist(c1, c2):
                unit_cycle_idx.remove(i)
                break
    final_cycles = [cycles[i] for i in unit_cycle_idx]
    result = []
    for item in final_cycles:
        result = result + item[:-1]
        # #print(result, "result")
    return result

def contains_sublist(A, B):
    if not B or not A:
        return False

    len_A, len_B = len(A), len(B)

    for i in range(len_A - len_B + 1):
        if A[i:i + len_B] == B:
            return True

    return False

def remove_cycle_in_final_all(twoDcycles, line_arr):
    cycle_arr = [[item.replace("+","").replace("-","") for item in cycle] for cycle in twoDcycles]
    line_arr = [item.replace("+","").replace("-","") for item in line_arr]
    line_set = set(line_arr)
    for item in cycle_arr:
        if set(item) == line_set:
            return True
    return False



# read the final.txt file, and filter
# kind 1: ABAB to AB
# kind 2: record A and record B have 90% similarity. keep the longest.
# kind 2:
def filter_final(final_all_file, cycle_count, cycle_result, ori_cycle_result, before_cut_dict):
    final_cycle_count = cycle_count
    final_all = open(final_all_file)
    line_idx = 0
    cycle_result_size = len(cycle_result)
    tmp_result = copy.deepcopy(cycle_result)
    before_cut_dict_swap = {v:k for k, v in before_cut_dict.items()}
    for line in final_all.readlines():
        if line.strip() == "":
            continue
        if line_idx < cycle_count:
            line_idx = line_idx + 1
            #continue
        line_k = line.strip().replace("\t","").replace("+","+\t").replace("-","-\t").strip()
        if line_k in before_cut_dict.keys():
            line_arr_tmp = re.split("\t", before_cut_dict[line_k])
        else:
            line_arr_tmp = re.split("\t", line_k)
        if remove_cycle_in_final_all(ori_cycle_result, line_arr_tmp):
            continue
        line_arr = re.split(r"\s+", line.strip())
        # for item in line_arr:
        #     sub_lst.append(item)
        # line_arr = line_arr.sort()
        # kind 1: for each record, filter cycle result like : ABABABAC to ABAC
        remove_duplicated_arr = line_arr_tmp 
        tmp_result.append(remove_duplicated_arr)
        line_idx = line_idx + 1
    keeped_idx = set(range(0, len(tmp_result)))
    for i in range(0, len(tmp_result)):
        if i not in keeped_idx:
            continue
        for j in range(i, len(tmp_result)):
            if i == j:
                continue
            if j not in keeped_idx:
                continue
            similar, idx = is_similar(tmp_result[i], tmp_result[j])
            if similar:
                if idx == 0:
                    # #print(j,"j")
                    keeped_idx.remove(j)
                    if j < cycle_count:
                        final_cycle_count=final_cycle_count-1
                else:
                    keeped_idx.remove(i)
                    if i < cycle_count:
                        final_cycle_count=final_cycle_count-1
                    break
    final_result = [tmp_result[i] for i in keeped_idx]
    final_result_cycle = []
    final_result_uncycle = []
    for item in final_result:
        if item in cycle_result:
            final_result_cycle.append(item)
        else:
            if "\t".join(item) in before_cut_dict_swap.keys():
                final_result_uncycle.append(before_cut_dict_swap["\t".join(item)].split("\t"))
            else:
                final_result_uncycle.append(item)
    return len(final_result_cycle), final_result_cycle+final_result_uncycle

def remove_dup_record(lst):
    lst_str = " ".join(lst)
    while True:
        out = re.sub(r'(?<!\S)(\S+(?:\s\S+)*)\s+\1(?!\S)', '\\1', lst_str)
        if out == lst_str:
            break
        lst_str = out
    return out.split(" ")
def remove_dup_record_shortest(a):
    n = len(a)
    final_list = a
    for i in range(n):
        rotated_lst = a[i:] + a[:i]
        current_list =remove_dup_record(rotated_lst)
        if len(current_list) < len(final_list):
            final_list = current_list
    return final_list



def make_fa_from_order(fain,orderin,faout,prefix,final_cycle_count):
    order_stream = open(orderin)
    faout_stream = open(faout, "w")
    record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
    n_seq = Seq("N" * 50)
    count = 0
    i = 0
    for line in order_stream.readlines():
        if "all" in line:
            continue
        seq = ""
        tmp = re.split("\t+", line.strip("\n"))
        for t in tmp:
            if t == '':
                continue
            t = t.replace("ref","")
            tmp_seq = record_dict[t[0:-1]].seq
            if t[-1] == '-':
                tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
            if seq == "":
                seq = seq  + tmp_seq
            else:
                seq = seq  + n_seq + tmp_seq
        count += 1
        tag = 'linear'
        if i < final_cycle_count:
            tag = 'cycle'
        i = i+1
        faout_stream.write(">" + prefix + "_phage_" + str(count)+"_"+tag + "\n")
        faout_stream.write(str(seq) + "\n")

    order_stream.close()
    faout_stream.close()

def get_seg_len(seg):
    return FAI_LEN[seg.replace("+", "").replace("-","")]

FAI_LEN = {}
if __name__ == "__main__":
    out_dir = sys.argv[1]
    prefix = sys.argv[2]
    cycle_file = sys.argv[3]
    final_all_file = sys.argv[4]
    filtered_final_all_file = os.path.join(out_dir, sys.argv[5])
    filtered_final_fa = os.path.join(out_dir, sys.argv[6])
    edge_fasta = sys.argv[7]
    cycle_out_file = os.path.join(out_dir, sys.argv[8])
    bam = sys.argv[9]
    before_cut = sys.argv[10]
    min_len = int(sys.argv[11])
    with open(edge_fasta+".fai", 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            FAI_LEN[fields[0]] = int(fields[1])
    os.makedirs(out_dir, exist_ok=True)
    
    before_cut_dict = {}
    with open(before_cut, "r") as f:
        for line in f:
            key,value=line.strip().split(":")
            if len(key) >0:
                before_cut_dict[key.strip()] = value.strip()
    # split edge fasta
    # prefix = get_prefix(Fasta(edge_fasta_file)) # SAMEA728574_phage_3, prefix:SAMEA728574
    # prefix = split_fasta(edge_fasta_file, final_all_file, edge_fasta)

    # remove_dup for cycle_file

    cycle_count, cycle_result, ori_cycle_result = filter_cycle(cycle_file, cycle_out_file,bam)

    # remove_dup for final.txt
    final_cycle_count, filtered_final_results = filter_final(final_all_file, cycle_count, cycle_result, ori_cycle_result, before_cut_dict)
    with open(filtered_final_all_file, "w") as filtered_final_all:
        for item in filtered_final_results:
            if get_path_len(item) > min_len:
                filtered_final_all.write("\t".join(item) + "\n")

    #make new fa
    make_fa_from_order(edge_fasta,filtered_final_all_file,filtered_final_fa,prefix,final_cycle_count)
