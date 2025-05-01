import sys
import re

fastg_fai = open(sys.argv[1])
original_graph = sys.argv[2]
outs = open(sys.argv[3], "w")
depth = int(float(sys.argv[4]))
f_th = sys.argv[5]
gene_file = sys.argv[6]
score_file = sys.argv[7]
blast_in = open(sys.argv[8])
blast_ratio = float(sys.argv[9])
fasta_fai = sys.argv[10]
hit_segs_out = open(sys.argv[11],"w")
contig_path = open(sys.argv[12])
to_remove_score_threhold = 0.2
relevate_edge_len = 200
fastg_graph = {}
gene_res = {}
scores = {}
sample = 'sss'
blast_segs = set()
relevate_blast_segs = set()
score_segs = set()
prev_seg = ""
prev_len = 0
prev_ref = ""
MIN_CYCLE_LEN = 1000
MIN_ALN_LEN = 1000

num_to_full_name = {}
full_name_to_num = {}
fai_len = {}
with open(fasta_fai, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        fai_len[fields[0]] = int(fields[1])
        f_0_list = fields[0].split("_")
        num_to_full_name[f_0_list[1]] = fields[0]
        full_name_to_num[fields[0]] = f_0_list[1]
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
for line in blast_in.readlines():
    t = line.strip().split("\t")
    #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
        #print(prev_len,"tdtd")
    if (prev_seg != t[0] and prev_seg != "") or (prev_ref != t[1] and prev_ref != ""):
        #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
            #print(prev_len, elen,"xdxd")
        elen = fai_len[prev_seg]
        if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)
        prev_seg = t[0]
        prev_ref = t[1]
        prev_len = int(t[3])
    else:
        if float(t[2]) > blast_ratio*100:
            #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
                #print(prev_len,float(t[2]),blast_ratio*100,"mmm")
            prev_len = prev_len + int(t[3])
        prev_seg = t[0]
        prev_ref = t[1]
if prev_seg != "":
    elen = fai_len[prev_seg]
    if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
        blast_segs.add(t[0])
# print(blast_segs)
if len(sys.argv) > 5:
    with open(gene_file, 'r') as gene_lst:
        for gene_r in gene_lst:
            contig_name,_ = gene_r.split("\t")
            gene_res[contig_name] = '1'

    with open(score_file, 'r') as score_lst:
        for score_r in score_lst:
            line_lst = score_r.strip().split("\t")
            scores[line_lst[0]] = str(line_lst[1])
            if float(line_lst[1]) > 0.7:
                score_segs.add(line_lst[0])
for line in fastg_fai:
    vs = line.split("\t")
    a = re.split(":|,|;", vs[0])
    fastg_graph[a[0]] = a[1:]
# print(tmp)
all_segs = {}
write_segs = set()
write_juncs = []
hit_segs = {}

with open(original_graph, 'r') as file:
    lines = file.readlines()
    
for line in lines:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        all_segs[vs[1]] = line
        # outs.write(line)
        if vs[1] not in hit_segs.keys():
            hit_segs[vs[1]] = ""
        if vs[1] in blast_segs:
                hit_segs[vs[1]] = hit_segs[vs[1]] + "ref+"
                relevate_blast_segs.add(vs[1])
                write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
        if float(scores[vs[1]]) > 0.7:
                hit_segs[vs[1]] = hit_segs[vs[1]] + "score+"
                relevate_blast_segs.add(vs[1])
                write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")                
        if vs[1] in gene_res:
                hit_segs[vs[1]] = hit_segs[vs[1]] + "gene+"
                write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
                relevate_blast_segs.add(vs[1])
        if hit_segs[vs[1]] != "":
            hit_segs_out.write(sample+"\t"+vs[1]+"\t"+hit_segs[vs[1]]+"\n")
        # outs.write(line)
        continue
    # if int(vs[-1]) < int(f_th):
    #     continue
    # print(vs)
    left_score = float(scores[vs[1]])
    right_score = float(scores[vs[3]])

    #if (left_score < 0.2 and left_node_len > 10000) or (right_score < 0.2 and right_node_len > 10000):
    #    continue
    k = vs[1]
    kc = vs[1] + "'"
    if vs[2] == '-':
        k = vs[1] + "'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3] + "'"
    if vs[4] == '-':
        v = vs[3] + "'"
        vc = vs[3]
    # if k =="EDGE_5369_length_3828_cov_7.082378":
    #     print(tmp[k])
    #     print(v)
    # print("k:", k)
    # print("kc", kc)
    # print("v", v)
    # print("vc", vc)
    is_add_junc = False

    if vs[1] == vs[3]:
        write_juncs.append(" ".join(str(item) for item in vs) + "\n")
        # print(all_segs[vs[1]].strip()+" "+(gene_res[vs[1]] if vs[1] in gene_res else '0')+" "+(scores[vs[1]] if vs[
        # 1] in scores else '0')+"\n")
        write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
        write_segs.add(all_segs[vs[3]].strip() + " " + (gene_res[vs[3]] if vs[3] in gene_res else '0') + " " + (
            str(scores[vs[3]]) if vs[3] in scores else '0') + " " + ('1' if vs[3] in blast_segs else '0') + "\n")
        relevate_blast_segs.add(vs[1])
        relevate_blast_segs.add(vs[3])

    if (vs[1] in blast_segs or vs[1] in gene_res or left_score >= 0.7) or (vs[3] in blast_segs or vs[3] in gene_res or right_score >= 0.7):

        # or (
            # vs[1] == vs[3]) or (left_score >= 0.6 and right_score >= 0.6):
        write_juncs.append(" ".join(str(item) for item in vs) + "\n")
        # print(all_segs[vs[1]].strip()+" "+(gene_res[vs[1]] if vs[1] in gene_res else '0')+" "+(scores[vs[1]] if vs[
        # 1] in scores else '0')+"\n")
        write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
        write_segs.add(all_segs[vs[3]].strip() + " " + (gene_res[vs[3]] if vs[3] in gene_res else '0') + " " + (
            str(scores[vs[3]]) if vs[3] in scores else '0') + " " + ('1' if vs[3] in blast_segs else '0') + "\n")
        relevate_blast_segs.add(vs[1])
        relevate_blast_segs.add(vs[3])
    else:
        if vs[1] in relevate_blast_segs or vs[3] in relevate_blast_segs:
            write_juncs.append(" ".join(str(item) for item in vs) + "\n")
            # print(all_segs[vs[1]].strip()+" "+(gene_res[vs[1]] if vs[1] in gene_res else '0')+" "+(scores[vs[1]] if vs[
            # 1] in scores else '0')+"\n")
            write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
                " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
            write_segs.add(all_segs[vs[3]].strip() + " " + (gene_res[vs[3]] if vs[3] in gene_res else '0') + " " + (
                str(scores[vs[3]]) if vs[3] in scores else '0') + " " + ('1' if vs[3] in blast_segs else '0') + "\n")
            relevate_blast_segs.add(vs[1])
            relevate_blast_segs.add(vs[3])
            


for line in lines:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        continue
    #if (left_score < 0.2 and left_node_len > 10000) or (right_score < 0.2 and right_node_len > 10000):
    #    continue
    k = vs[1]
    kc = vs[1] + "'"
    if vs[2] == '-':
        k = vs[1] + "'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3] + "'"
    if vs[4] == '-':
        v = vs[3] + "'"
        vc = vs[3]
    if vs[1] in relevate_blast_segs or vs[3] in relevate_blast_segs:
        write_juncs.append(" ".join(str(item) for item in vs) + "\n")
        # print(all_segs[vs[1]].strip()+" "+(gene_res[vs[1]] if vs[1] in gene_res else '0')+" "+(scores[vs[1]] if vs[
        # 1] in scores else '0')+"\n")
        write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + str(scores[vs[1]]) if vs[1] in scores else '0') + " " + ('1' if vs[1] in blast_segs else '0') + "\n")
        write_segs.add(all_segs[vs[3]].strip() + " " + (gene_res[vs[3]] if vs[3] in gene_res else '0') + " " + (
            str(scores[vs[3]]) if vs[3] in scores else '0') + " " + ('1' if vs[3] in blast_segs else '0') + "\n")


blast_segs.update(set(gene_res.keys()))
blast_segs.update(score_segs)
support_segs = blast_segs
path_segs = filter_paths(contig_path, num_to_full_name, support_segs)
# print("xxxxx",len(gene_res),len(score_segs),len(blast_segs),len(write_segs),len(path_segs))
pure_segs = []
for item in write_segs:
    pure_seg = item.split(" ")[1]
    pure_segs.append(pure_seg)
    # print(item)
    outs.write(item)
for seg in path_segs:
    # for seg in path:
    if seg not in pure_segs:
        outs.write(f"{all_segs[seg].strip()} 0 1.0 0\n")
seen_juncs = set()
for item in write_juncs:
    if item not in seen_juncs:
        outs.write(item)
        seen_juncs.add(item)
