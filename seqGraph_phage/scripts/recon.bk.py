import sys
from operator import itemgetter
import re
import subprocess

def run_samtools_depth(input_bam, region):
    cmd = f'samtools depth -r {region} {input_bam}'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Error running samtools depth: {result.stderr}")
        return None

    depths = [int(line.split('\t')[2]) for line in result.stdout.splitlines()]
    return depths

def calculate_average_depth(depths):
    total_depth = sum(depths)
    average_depth = total_depth / len(depths)
    return average_depth
def get_depth(segs, bam):
    total_depths = []
    total_positions = 0
    average_contig_depth = {}
    final_segs = []
    for item in segs:
        contig = list(item.split())[1]
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
        i_arr[-1] = str(copy_num)
        final_segs.append(" ".join(i_arr))
    return final_segs

#fai = open(sys.argv[1])  # fai
graph = open(sys.argv[1])  # graph
segf = open(sys.argv[3])  # 需要从新的segs
outs = sys.argv[2]  # 输出prefix
depth = float(sys.argv[4])  # samtools 计算而来
min_support=1
if len(sys.argv) > 5:
    min_support = int(sys.argv[5])
if len(sys.argv) > 6:
    bam = sys.argv[6]
tmp = {}
segs = []
seg_g = {}
out_juncs = []
#for line in fai:
#    vs = line.split("\t")
#    a = re.split(":|,|;", vs[0])
#    tmp[a[0]] = a[1:]
ref_name_records = {}
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
for js in out_segs:
    outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "w")
    i = i + 1
    js_c = get_depth(js,bam)
    # print(i, js)
    for j in js_c:
        outs_f.write(j+"\n")
    outs_f.close()
i = 0
# print(out_juncs)
for js in out_juncs:
    outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "a")
    i = i + 1
    for j in sorted(js):
        if int(j[-1]) < min_support:
            continue
        outs_f.write(" ".join(j) + "\n")
    outs_f.close()
