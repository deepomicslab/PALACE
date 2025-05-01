import re
import os
import sys
res = set()
ignore_len = int(sys.argv[2])
gene_hit = sys.argv[3]
score_file = sys.argv[4]
min_gene_count_for_single_contig = 10
with open(sys.argv[1]) as r:
    for line in r.readlines():
        if "loop" in line or "iter" in line:
            continue
        line_len = 0
        splited=re.split(r'[+-]',line.strip())
        for v in splited:
            if v == "" or v ==" ":
                continue
            if ignore_len == 0:
                line_len = line_len + int(v.split('_')[3])
        if line_len >= 10000 and ignore_len == 0:
        #print(line,"xxx\n")
            liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","").replace("ref","")
            res.add(liner.strip("\n"))
        else:
            liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","").replace("ref","")
            res.add(liner.strip("\n"))
genehit = []
scorehit = []
with open(gene_hit, 'r') as gh:
    for s in gh:
        item = s.strip().split('\t')
        if len(s.strip()) == 0:
            continue
        if int(item[1]) >= min_gene_count_for_single_contig:
            genehit.append(item[0])
with open(score_file, 'r') as ps:
    for s in ps:
        item = s.strip().split('\t')
        if float(item[1]) >= 0.7:
            scorehit.append(item[0])
for item in res:
    item.replace("+","+\t").replace("-","-\t")
    result_item = item
    item_list = item.split("\t")
    if len(item_list) <= 1:
        if item_list[0] not in genehit and item_list[0] not in scorehit:
            continue
    print(result_item)
