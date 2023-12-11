import sys
from operator import itemgetter
import re
import subprocess


def parse_blast(blast_in, blast_ratio, refs):
    ref_blast_segs = {}
    prev_seg = ""
    prev_len = 0
    prev_ref = ""
    prev_seg_len = 0
    for ref in refs:
        ref_blast_segs[ref] = set()
    for line in blast_in.readlines():
        t = line.strip().split("\t")
        #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
            #print(prev_len,"tdtd")
        if t[1] not in refs:
            continue
        if (prev_seg != t[0] and prev_seg != "") or (prev_ref != t[1] and prev_ref != ""):
            #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
                #print(prev_len, elen,"xdxd")
            elen = prev_seg_len
            if float(prev_len) / float(elen) > blast_ratio:
                ref_blast_segs[prev_ref].add(prev_seg)
                #blast_segs.add(prev_seg)
            prev_seg = t[0]
            prev_ref = t[1]
            prev_len = int(t[5])
            prev_seg_len = int(t[3])
        else:
            if float(t[2]) > blast_ratio*100:
                #if "EDGE_27711765_length_1685_cov_5" in prev_seg:
                    #print(prev_len,float(t[2]),blast_ratio*100,"mmm")
                prev_len = prev_len + int(t[5])
            prev_seg = t[0]
            prev_ref = t[1]
            prev_seg_len = int(t[3])
    if prev_seg != "":
        elen = prev_seg_len
        if float(prev_len) / float(elen) > blast_ratio:
            ref_blast_segs[prev_ref].add(prev_seg)
            #blast_segs.add(t[0])
    return ref_blast_segs


def run_samtools_depth(input_bam, region):
    cmd = f'{samtools} depth -r {region} {input_bam}'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Error running samtools depth: {result.stderr}")
        return None

    depths = [int(line.split('\t')[2]) for line in result.stdout.splitlines()]
    return depths

def calculate_average_depth(depths):
    total_depth = sum(depths)
    if len(depths) == 0:
        return 0
    average_depth = total_depth / len(depths)
    return average_depth
def get_depth(segs, bam):
    total_depths = []
    total_positions = 0
    average_contig_depth = {}
    final_segs = []
    for item in segs:
        contig = item
        depths = run_samtools_depth(bam, contig)
        if depths:
            average_depth = calculate_average_depth(depths)
            average_contig_depth[item] = average_depth
            total_depths.extend(depths)
            total_positions += len(depths)
    total_average_depth = calculate_average_depth(total_depths)
    for k in average_contig_depth.keys():
        i_arr = k.split()
        i_depth = average_contig_depth[k]
        copy_num = round(i_depth/total_average_depth)
        if copy_num == 0:
            copy_num = 1
        final_segs.append("SEG " + k + " " +k.split('_')[-1] + " "+str(copy_num))
    return final_segs

#fai = open(sys.argv[1])  # fai
graph = open(sys.argv[1])  # graph
segf = open(sys.argv[3])  # 需要从新的segs
outs = sys.argv[2]  # 输出prefix
samtools = sys.argv[4]  # samtools
min_support=1
if len(sys.argv) > 5:
    min_support = int(sys.argv[5])
if len(sys.argv) > 6:
    bam = sys.argv[6]
if len(sys.argv) > 7:
    blast_ratio = float(sys.argv[7])
tmp = {}
segs = []
seg_g = {}
out_juncs = []
#for line in fai:
#    vs = line.split("\t")
#    a = re.split(":|,|;", vs[0])
#    tmp[a[0]] = a[1:]
ref_name_records = {}
ref_segs = {}
for idx, line in enumerate(segf):
    l = line.strip().split('\t')
    t = re.split(r"[\+\-]", l[0])[:-1]
    ref_name_records[idx] = l[1]
    line_segs = []
    out_juncs.append([])
    for v in t:
        line_segs.append(v)
        if v in seg_g.keys():
            if idx in seg_g[v]:
                continue
            seg_g[v].append(idx)
        else:
            seg_g[v] = [idx]
    segs.append(line_segs)
    ref_segs[l[1]] = set(t)
out_segs = [[] for i in range(0, len(out_juncs))]


def in_segs(seg1, seg2=None):
    if seg2 is None:
        for s in segs:
            if seg1 in s:
                return 1
    else:
        res = []
        for i in range(len(segs)):
            s = segs[i]
            if seg1 in s and seg2 in s:
                res.append(i+1)
        return res
    return 0


for line in graph:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        if in_segs(vs[1]):
            # print(vs, 'test', vs[1])
            for v in seg_g[vs[1]]:
                out_segs[v].append(line)
        continue

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
    # if vs[1] == "EDGE_72432_length_111_cov_3.857143" or vs[3] =="EDGE_72432_length_111_cov_3.857143":
        # print(in_segs(vs[1], vs[3]), "rrrrrrrrrr")
    # if vs ==[]:
        # print(in_segs(vs[1], vs[3]), 'dddddddd')

    if in_segs(vs[1], vs[3]):
        for it in in_segs(vs[1], vs[3]):
            out_juncs[it - 1].append(vs)
# for item in tmp.keys():
#     first = item
#     fdir = "+"
#     if item[-1] == "'":
#         first = item[0:-1]
#         fdir = "-"
#     for i in tmp[item]:
#         if i == "":
#             continue
#         second = i
#         sdir = "+"
#         if i[-1] == "'":
#             second = i[0:-1]
#             sdir = "-"
#         if in_segs(first, second):
#             vs = ["JUNC",first, fdir, second, sdir, str(int(depth) / 2)]
#             out_juncs[in_segs(vs[1], vs[3]) - 1].append(vs)
            # out_juncs.append(vs)
            # out_juncs("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
# s_juncs = sorted(out_juncs, key=itemgetter(2))
i = 0
# print(out_segs)
# filter segs mapped to reof
#blast_stream = open(blast_in)
#ref_segs = parse_blast(blast_stream, blast_ratio, list(ref_name_records.values()))
#for js in out_segs:
    #outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "w")
    #js_c = get_depth(js,bam)
#    ok_segs = parse_blast(blast_stream, blast_ratio, ref_name_records[i])
#    print(ok_segs)
#    ref_segs[ref_name_records[i]] = set(ok_segs)
#    i = i + 1
    # print(i, js)
    #for j in j_c:
    #   if j.split()[1] in ok_segs:
    #        outs_f.write(j+"\n")
    #outs_f.close()
i = 0
# print(out_juncs)
for idx,js in enumerate(out_juncs):
    outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "w")
    current_count = len(ref_segs[ref_name_records[i]])
    prev_ref_segs_count = 0
    while prev_ref_segs_count != current_count:
        prev_ref_segs_count = current_count
        for j in sorted(js):
            if int(j[-1]) < min_support or (j[1] not in ref_segs[ref_name_records[i]] and j[3] not in ref_segs[ref_name_records[i]]):
                continue
            ref_segs[ref_name_records[i]].add(j[1])
            ref_segs[ref_name_records[i]].add(j[3])
        current_count = len(ref_segs[ref_name_records[i]])
    if len(js) == 0:
        ref_segs[ref_name_records[i]] = segs[idx]
    ref_segs_copy_corrected = get_depth(ref_segs[ref_name_records[i]], bam)
    for seg in ref_segs_copy_corrected:
        outs_f.write(seg+"\n")
    for j in sorted(js):
        if int(j[-1]) < min_support or (j[1] not in ref_segs[ref_name_records[i]] and j[3] not in ref_segs[ref_name_records[i]]):
            continue
        outs_f.write(" ".join(j)+"\n")
    i = i + 1
    outs_f.close()