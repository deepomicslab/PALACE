import sys
from pyfaidx import Fasta
import re
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
import copy
import numpy as np
from collections import Counter, defaultdict
from itertools import chain

# ============================================================
# Global cache for samtools depth results to avoid redundant calls
# ============================================================
SAMTOOLS_RAW_DEPTH_CACHE = {}

def get_cached_raw_depths(contig, bam):
    """
    Cached samtools depth query returning raw depth array.
    Returns cached result if available, otherwise calls samtools and caches.
    """
    if contig in SAMTOOLS_RAW_DEPTH_CACHE:
        return SAMTOOLS_RAW_DEPTH_CACHE[contig]
    depths = run_samtools_depth(bam, contig)
    SAMTOOLS_RAW_DEPTH_CACHE[contig] = depths  # May be None
    return depths

# ============================================================
# Parser and smart quota deduplication (using node cov values)
# ============================================================

def parse_line_nodes(line):
    """
    Parse node information from a line, extracting coverage values from node names.
    """
    pattern = re.compile(r'(EDGE_(\d+)_length_(\d+)_cov_([\d\.]+)([+-]))')
    matches = pattern.findall(line)

    nodes = []
    for m in matches:
        try:
            nodes.append({
                'full': m[0],       # Full node string
                'id': m[1],         # Node ID
                'len': int(m[2]),   # Node length
                'cov': float(m[3])  # Coverage value from node name
            })
        except ValueError as e:
            print(f"Warning: Could not parse node: {m[0]}. Error: {e}")
            continue

    return nodes

def calculate_baseline(nodes):
    """
    Calculate baseline single-copy coverage.
    Preferentially uses median coverage of nodes appearing once in the path.
    """
    if not nodes:
        return 1.0

    id_counts = Counter([n['id'] for n in nodes])
    single_copy_covs = [n['cov'] for n in nodes if id_counts[n['id']] == 1]

    if single_copy_covs:
        return np.median(single_copy_covs)
    else:
        return np.median([n['cov'] for n in nodes])

def smart_quota_dedup(line, bam=None):
    """
    Coverage quota-based deduplication logic using node coverage values.
    """
    line = line.strip()
    if not line:
        return ""

    nodes = parse_line_nodes(line)
    if not nodes:
        return line

    baseline = calculate_baseline(nodes)
    if baseline == 0:
        baseline = 1.0

    is_hub = lambda cov: cov > 2.5 * baseline

    cov_by_id = {}
    for n in nodes:
        uid = n['id']
        cov_by_id[uid] = max(cov_by_id.get(uid, 0.0), n['cov'])

    node_budget = defaultdict(int)
    for uid, max_cov in cov_by_id.items():
        if is_hub(max_cov):
            node_budget[uid] = 999999
        else:
            copies = int(round(max_cov / baseline))
            node_budget[uid] = max(1, copies)

    temp_path = []
    for node in nodes:
        uid = node['id']
        if node_budget[uid] > 0:
            temp_path.append(node)
            node_budget[uid] -= 1

    if not temp_path:
        return ""

    final_path_strings = []
    last_node_full = None
    for node in temp_path:
        current_full = node['full']
        if current_full != last_node_full:
            final_path_strings.append(current_full)
            last_node_full = current_full

    return "\t".join(final_path_strings)

def apply_smart_quota_dedup_to_path(path_list, bam):
    """
    Apply smart_quota_dedup to a path list.
    Returns deduplicated path list.
    """
    line_str = "\t".join(path_list)
    deduped_str = smart_quota_dedup(line_str, bam)
    if not deduped_str:
        return []
    return deduped_str.split("\t")


# ============================================================
# Original corrected_dup.py functions (cycle logic uses samtools depth)
# ============================================================

def get_path_len(path):
    sum = 0
    for item in path:
        if item.startswith("EDGE"):
            arr = item.split("_")
            sum += int(arr[3])
    return sum

def split_list_on_element(input_list, A):
    indices = [i for i, elem in enumerate(input_list) if A in elem]
    indices.append(len(input_list))
    sublists = [input_list[indices[i]:indices[i + 1]] for i in range(len(indices) - 1)]
    sublist_counts = Counter(tuple(sublist) for sublist in sublists)
    return sublist_counts

def concatenate_sublists_with_counts(sublist_counts):
    repeated_sublists = [list(sublist) * count for sublist, count in sublist_counts.items()]
    flattened_list = list(chain.from_iterable(repeated_sublists))
    return flattened_list

def merge_repeat(lst):
    s = [item.replace("-", "").replace("+", "") for item in lst]
    item_counts = Counter(s)
    most_frequent_item = max(item_counts, key=item_counts.get)
    item_index = s.index(most_frequent_item)
    reformatted_list = lst[item_index:] + lst[:item_index]
    sublist_counts = split_list_on_element(reformatted_list, most_frequent_item)
    return concatenate_sublists_with_counts(sublist_counts)

def run_samtools_depth(input_bam, region):
    cmd = f'samtools depth -r {region} {input_bam}'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
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
    min_seg = ""
    min_copy = 10000
    for item in unit_seg:
        item = item.replace("+", "").replace("-", "")
        if item not in seg_copies.keys():
            copy = 1
        else:
            copy = seg_copies[item]
        if copy < min_copy:
            min_seg = item
            min_copy = copy
    return min_seg, min_copy

def non_dup_item(ori_arr, unit_cycles):
    ori_arr_str = "\t".join(ori_arr).replace("+", "").replace("-", "")
    unit_cycle_str = ["\t".join(item).replace("+", "").replace("-", "") for item in unit_cycles]
    for item in unit_cycle_str:
        ori_arr_str.replace(item, "")
    return ori_arr_str.split("\t")

def calculate_real_copy_for_cycle(unit_seg, seg_copies, non_unit_part):
    min_seg, min_copy = get_min_copy_seg(unit_seg, seg_copies)
    other_count = non_unit_part.count(min_seg)
    real_copy = min_copy - other_count
    if real_copy < 1:
        real_copy = 1
    return real_copy

def get_depth(all_non_dup_segs, unit_cycle_segs, non_unit_part, bam, first_item):
    """
    Calculate depth for cycle-related logic using cached samtools depth.
    """
    unit_copies = []
    seg_len_depth = {}
    total_depths = []
    total_positions = 0
    seg_depth_dict = {}

    for item in all_non_dup_segs:
        contig = item.replace('-', '').replace('+', '')
        depths = get_cached_raw_depths(contig, bam)
        if depths:
            average_depth = calculate_average_depth(depths)
            seg_len_depth[contig] = [average_depth, len(depths)]
            total_depths.extend(depths)
            total_positions += len(depths)

    total_average_depth = calculate_average_depth(total_depths)
    for k in seg_len_depth.keys():
        if total_average_depth > 0:
            seg_depth_dict[k] = round(seg_len_depth[k][0] / total_average_depth)
        else:
            seg_depth_dict[k] = 1

    for unit_seg in unit_cycle_segs:
        seg_unit_copy = calculate_real_copy_for_cycle(unit_seg, seg_depth_dict, non_unit_part)
        unit_copy = round(seg_unit_copy)
        if unit_copy == 0:
            unit_copy = 1
        unit_copies.append(unit_copy)

    k = first_item.replace('-', '').replace('+', '')
    return_depth = 0
    if k in seg_depth_dict.keys():
        return_depth = seg_depth_dict[k]
    return unit_copies, return_depth

def reformat_cycle(s):
    ori_s = copy.deepcopy(s)
    n = len(s)
    longest_substring_index = -1
    for i in range(n // 2 + 1):
        if s[:i] == s[-i:]:
            longest_substring_index = i
    if longest_substring_index != -1:
        return s[len(s) - longest_substring_index:] + s[0: len(s) - longest_substring_index]
    if ori_s == s:
        s = merge_repeat(ori_s)
    return s

def are_cyclically_equal(s1, s2):
    if s1 in s2:
        return True
    concatenated = s1 + "\t" + s1
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
                    if are_cyclically_equal(item, "\t".join(s[start:start + repeat_len])):
                        tag = False
                        break
                if tag:
                    repeats.add("\t".join(s[start:start + repeat_len]))
    return [item.split("\t") for item in repeats]

def get_prefix(edge_fasta):
    return list(edge_fasta.keys())[0].split("_")[0]

def split_fasta(edge_fasta_file, final_all_file, edge_fasta_out_file):
    edge_fasta = Fasta(edge_fasta_file)
    edge_txt = open(final_all_file)
    edge_out = open(edge_fasta_out_file, "w")
    writed_edge = []
    idx = 1
    prefix = get_prefix(edge_fasta)
    for line in edge_txt.readlines():
        prev_len = 0
        line_arr = re.split(r'\s+', line.strip())
        edge_fasta_key = prefix + "_phage_" + str(idx)
        for item in line_arr:
            item_dir = item[-1]
            item_len = get_seg_len(item)
            item_fa = edge_fasta[edge_fasta_key][prev_len:item_len + prev_len]
            if item_dir == "-":
                item_fa_rev = item_fa.reverse.complement
            prev_len = item_len + prev_len
            if item[:-1] in writed_edge:
                continue
            edge_out.write(">" + item[:-1] + "\n")
            writed_edge.append(item[:-1])
            if item_dir == "-":
                edge_out.write(str(item_fa_rev) + "\n")
            else:
                edge_out.write(str(item_fa) + "\n")
        idx = idx + 1
    edge_txt.close()
    edge_out.close()
    return prefix

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

def count_item_ignor_direction(lst, ele):
    count = 0
    ele = ele.replace("+", "").replace("-", "")
    for item in lst:
        if ele in item:
            count = count + 1
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
        if unit_copy < 1:
            unit_copy = 1
        start_idx, end_idx = find_sublist_indexes(unit_item, line_arr)
        line_arr = line_arr[:start_idx] + unit_cycles[i] * unit_copy + line_arr[end_idx:]
    first_item_count_in_liner_arr = count_item_ignor_direction(line_arr, line_arr[0])
    if abs(first_item_count_in_liner_arr - first_item_copy) <= 1:
        return line_arr
    sublist_counts = split_list_on_element(line_arr, line_arr[0])
    final_list = []
    final_list_len = 0
    for sublist, count in sublist_counts.items():
        current_len = get_contig_len_for_arr(sublist)
        if current_len > final_list_len:
            final_list = sublist
            final_list_len = current_len
    return final_list

def filter_cycle(cycle_file, cycle_out_file, bam):
    cycle_out = open(cycle_out_file, "w")
    tmp_result = []
    line_count = 0
    ori_cycle_result = []
    for line in open(cycle_file).readlines():
        line_count = line_count + 1
        line_arr = re.split(r"\s+", line.strip())
        ori_cycle_result.append(line_arr)
        line_arr = reformat_cycle(line_arr)
        first_item = line_arr[0]
        unit_cycles = find_consecutive_repeats(line_arr)
        non_unit_part = non_dup_item(line_arr, unit_cycles)
        unit_copies, first_item_copy = get_depth(set(line_arr), unit_cycles, non_unit_part, bam, first_item)
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
    # If needed, uncomment to write cycle_out_file
    # for row in final_result:
    #     line = ' '.join(str(elem) for elem in row) + '\n'
    #     cycle_out.write(line)
    # cycle_out.close()
    return line_count, final_result, ori_cycle_result

def is_same(lst1, lst2):
    return set(lst2).issubset(set(lst1))

def is_similar(lst1, lst2):
    lst1_lens = [get_seg_len(item) for item in lst1]
    lst2_lens = [get_seg_len(item) for item in lst2]
    lst1_lens_sum = sum(set(lst1_lens))
    lst2_lens_sum = sum(set(lst2_lens))
    intersections = set(lst1_lens).intersection(lst2_lens)
    if sum(intersections) / lst1_lens_sum >= 0.9 or sum(intersections) / lst2_lens_sum >= 0.9:
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
    for i in range(0, len(cycles)):
        c1 = cycles[i]
        for j in range(i, len(cycles)):
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
    return result

def contains_sublist(A, B):
    if not B or not A:
        return False
    len_A, len_B = len(A), len(B)
    for i in range(len_A - len_B + 1):
        if A[i:i + len(B)] == B:
            return True
    return False

def remove_cycle_in_final_all(twoDcycles, line_arr):
    cycle_arr = [[item.replace("+", "").replace("-", "") for item in cycle] for cycle in twoDcycles]
    line_arr = [item.replace("+", "").replace("-", "") for item in line_arr]
    line_set = set(line_arr)
    for item in cycle_arr:
        if set(item) == line_set:
            return True
    return False

def filter_final(final_all_file, cycle_count, cycle_result, ori_cycle_result, before_cut_dict):
    final_cycle_count = cycle_count
    final_all = open(final_all_file)
    line_idx = 0
    cycle_result_size = len(cycle_result)
    tmp_result = copy.deepcopy(cycle_result)
    before_cut_dict_swap = {v: k for k, v in before_cut_dict.items()}
    for line in final_all.readlines():
        if line.strip() == "":
            continue
        if line_idx < cycle_count:
            line_idx = line_idx + 1
        line_k = line.strip().replace("\t", "").replace("+", "+\t").replace("-", "-\t").strip()
        if line_k in before_cut_dict.keys():
            line_arr_tmp = re.split("\t", before_cut_dict[line_k])
        else:
            line_arr_tmp = re.split("\t", line_k)
        if remove_cycle_in_final_all(ori_cycle_result, line_arr_tmp):
            continue
        line_arr = re.split(r"\s+", line.strip())
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
                    keeped_idx.remove(j)
                    if j < cycle_count:
                        final_cycle_count = final_cycle_count - 1
                else:
                    keeped_idx.remove(i)
                    if i < cycle_count:
                        final_cycle_count = final_cycle_count - 1
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
    return len(final_result_cycle), final_result_cycle + final_result_uncycle

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
        current_list = remove_dup_record(rotated_lst)
        if len(current_list) < len(final_list):
            final_list = current_list
    return final_list

def make_fa_from_order(fain, orderin, faout, prefix, final_cycle_count):
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
            t = t.replace("ref", "")
            tmp_seq = record_dict[t[0:-1]].seq
            if t[-1] == '-':
                tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
            if seq == "":
                seq = seq + tmp_seq
            else:
                seq = seq + n_seq + tmp_seq
        count += 1
        tag = 'linear'
        if i < final_cycle_count:
            tag = 'cycle'
        i = i + 1
        faout_stream.write(">" + prefix + "_phage_" + str(count) + "_" + tag + "\n")
        faout_stream.write(str(seq) + "\n")

    order_stream.close()
    faout_stream.close()

def get_seg_len(seg):
    return FAI_LEN[seg.replace("+", "").replace("-", "")]


# ============================================================
# Main program
# ============================================================

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

    with open(edge_fasta + ".fai", 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            FAI_LEN[fields[0]] = int(fields[1])

    os.makedirs(out_dir, exist_ok=True)

    before_cut_dict = {}
    with open(before_cut, "r") as f:
        for line in f:
            key, value = line.strip().split(":")
            if len(key) > 0:
                before_cut_dict[key.strip()] = value.strip()

    # Step 1: Process cycle file
    cycle_count, cycle_result, ori_cycle_result = filter_cycle(cycle_file, cycle_out_file, bam)

    # Step 2: Process final.txt file with deduplication and similarity filtering
    final_cycle_count, filtered_final_results = filter_final(
        final_all_file, cycle_count, cycle_result, ori_cycle_result, before_cut_dict
    )

    # Step 3: Apply smart_quota_dedup to each path
    quota_deduped_results = []
    for idx, path in enumerate(filtered_final_results):
        deduped_path = apply_smart_quota_dedup_to_path(path, bam)
        if deduped_path:
            quota_deduped_results.append(deduped_path)
        else:
            quota_deduped_results.append(path)

    # Step 4: Write final results to file
    with open(filtered_final_all_file, "w") as filtered_final_all:
        for item in quota_deduped_results:
            if get_path_len(item) > min_len:
                filtered_final_all.write("\t".join(item) + "\n")