import sys
import re

original_result = open(sys.argv[1])
outs = open(sys.argv[2], "w")
blast_in = open(sys.argv[3])
fasta_fai = sys.argv[4]
cycles_file = sys.argv[5]
blast_segs = set()
fai_len = {}
scores = {}
gene_res = []
def get_len(edge):
    return fai_len[edge]
def get_edge_len(edge):
    return int(edge.split("_")[3])
def filter_paths(contig_paths, num_to_full_name, support_segs):
    full_name_results = set()
    index = 1
    for line in contig_paths.readlines():
        line = line.strip().replace(";","")
        if line.startswith("NODE"):
            continue
        # if index % 2 == 0 and index % 4 != 0:
        full_names = []
        is_add = False
        nums = line.split(",")
        full_len = 0
        add_len = 0
        for num in nums:
            full_name = num_to_full_name[num[:-1]]
            full_names.append(full_name)
            e_len = get_edge_len(full_name)
            full_len = full_len + e_len
            if full_name in support_segs:
                is_add = True
                add_len = add_len + e_len
        if is_add:
            if add_len/full_len >= 0.5 or add_len > 2000:
                for n in full_names:
                    full_name_results.add(n)
    index = index + 1
    return full_name_results
cycle_segs = []
with open(cycles_file) as cf:
    for line in cf.readlines():
        line_replace = line.replace("+","").replace("-","")
        if "iter" in line_replace or "self" in line_replace:
            continue
        line_arr = line_replace.strip().split()
        cycle_segs.extend(line_arr)
prev_seg = ""
prev_ref = ""
prev_len = 0
fai_len = {}
with open(fasta_fai, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        fai_len[fields[0]] = int(fields[1])
for line in blast_in.readlines():
    t = line.strip().split("\t")
    #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
        #print(prev_len,"tdtd")
    if (prev_seg != t[0] and prev_seg != "") or (prev_ref != t[1] and prev_ref != ""):
        #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
            #print(prev_len, elen,"xdxd")
        elen = fai_len[prev_seg]
        if float(prev_len) / float(elen) > 0.7 or prev_len > 2000:
            blast_segs.add(prev_seg)
        prev_seg = t[0]
        prev_ref = t[1]
        prev_len = int(t[3])
    else:
        if float(t[2]) > 0.7*100:
            #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
                #print(prev_len,float(t[2]),blast_ratio*100,"mmm")
            prev_len = prev_len + int(t[3])
        prev_seg = t[0]
        prev_ref = t[1]
if prev_seg != "":
    elen = fai_len[prev_seg]
    if float(prev_len) / float(elen) > 0.7 or prev_len > 2000:
        blast_segs.add(t[0])
print(blast_segs)
# print(tmp)
all_segs = {}
write_segs = set()
write_juncs = []
hit_segs = {}
for line in original_result:
    ref_len = 0
    total_len = 0
    if "iter" in line or "self" in line:
        continue
    line_replaced = line.rstrip().replace("+","").replace("-","")
    split_text = re.split(r'\t+', line_replaced)
    # print(vs)
    is_cycle = False
    for item in split_text:
        if item in cycle_segs:
            is_cycle = True
            break
        item_len = get_edge_len(item)
        total_len = total_len + item_len
        if item in blast_segs:
            ref_len = ref_len + item_len
            # is_ok = True
            # break
    is_ok = False
    print(ref_len, total_len, is_cycle)
    if not is_cycle and ref_len/total_len >= 0.7:
        is_ok = True
    if is_ok:
        outs.write(line)
