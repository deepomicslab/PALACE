import argparse
import matplotlib.pyplot as plt
from operator import itemgetter
import matplotlib.patches as patches
import re
from collections import defaultdict
import os
def determine_strand_for_pair(file_path, query_of_interest, reference_of_interest):
    strand_lengths = defaultdict(int)
    with open(file_path, 'r') as file:
        for line in file:
            tokens = line.split()
            query, reference, qstart, qend, sstart, send = tokens[0], tokens[1], int(tokens[8]), int(tokens[9]), int(tokens[10]), int(tokens[11])
            if query == query_of_interest and reference == reference_of_interest:
                alignment_length = abs(qend - qstart) + 1
                if sstart < send:
                    strand_lengths['+'] += alignment_length
                else:
                    strand_lengths['-'] += alignment_length

    if strand_lengths['+'] > strand_lengths['-']:
        return '+'
    else:
        return '-'
def conver_minus_strand2_plus(query_name,cut_pos,fai_len):
    start_query = re.split(r'(\+|-)',query_name)
    start_query = [start_query[n] + start_query[n + 1] for n in range(0, len(start_query) - 1, 2)]
    total_len = get_line_len(query_name, fai_len)
    results_query = ""
    for item in reversed(start_query):
        if item[-1] == "-":
            new_item = item[:-1] + "+"
        else:
            new_item = item[:-1] + "-"
        results_query += new_item
    return results_query,total_len - cut_pos
def cut_end_contig(input_blast, blast_segs, fai_len,ref):
    ref_query_dict = defaultdict(lambda: {
        "min_start": float('inf'), "min_start_query": "",
        "max_end": float('-inf'), "max_end_query": "",
        "min_start_query_start": 0, "min_start_query_end": 0,
        "max_end_query_start": 0, "max_end_query_end": 0
    })

    with open(input_blast, "r") as file:
        lines = file.readlines()
        for line in lines:
            parts = line.strip().split("\t")
            query = parts[0]
            if query not in blast_segs:
                continue
            if query not in blast_segs:
                continue
            reference = parts[1]
            if reference not in ref:
                continue
            sstart = min(int(parts[10]),int(parts[11]))
            send = max(int(parts[11]),int(parts[10]))
            qstart = min(int(parts[8]),int(parts[9]))
            qend = max(int(parts[9]),int(parts[8]))

            #print(sstart, ref_query_dict[reference]["min_start"])
            if sstart < ref_query_dict[reference]["min_start"] or ref_query_dict[reference]["min_start_query"] == query:
                if ref_query_dict[reference]["min_start_query"] != query:
                    ref_query_dict[reference]["min_start"] = sstart
                    ref_query_dict[reference]["min_start_query"] = query
                    ref_query_dict[reference]["min_start_query_start"] = qstart
                    ref_query_dict[reference]["min_start_query_end"] = qend
                else:
                    ref_query_dict[reference]["min_start"] = sstart
                    if ref_query_dict[reference]["min_start_query_start"] > qstart:
                        ref_query_dict[reference]["min_start_query_start"] = qstart
                    if ref_query_dict[reference]["min_start_query_end"] < qend:
                        ref_query_dict[reference]["min_start_query_end"] = qend

            if send > ref_query_dict[reference]["max_end"] or ref_query_dict[reference]["max_end_query"] == query:
                if ref_query_dict[reference]["max_end_query"] != query:
                    ref_query_dict[reference]["max_end"] = send
                    ref_query_dict[reference]["max_end_query"] = query
                    ref_query_dict[reference]["max_end_query_start"] = qstart
                    ref_query_dict[reference]["max_end_query_end"] = qend
                else:
                    ref_query_dict[reference]["max_end"] = send
                    if ref_query_dict[reference]["max_end_query_end"] < qend:
                        ref_query_dict[reference]["max_end_query_end"] = qend
                    if ref_query_dict[reference]["max_end_query_start"] > qstart:
                        ref_query_dict[reference]["max_end_query_start"] = qstart
    # For each reference print the query split by start and end
    # cut by blast
    ref_start_end_segs={}
    for ref, data in ref_query_dict.items():
        strand = determine_strand_for_pair(input_blast, data["min_start_query"],ref)
        original_min_start_query = data["min_start_query"]
        if strand == "-":
            data["min_start_query"],data["min_start_query_start"] = conver_minus_strand2_plus(data["min_start_query"], data["min_start_query_end"],fai_len)
        #print("start",strand, data["min_start_query"], data["min_start_query_start"])
        start_query = re.split(r'(\+|-)',data["min_start_query"])
        start_query = [start_query[n] + start_query[n + 1] for n in range(0, len(start_query) - 1, 2)]
        start_query_start = data["min_start_query_start"]
        #start_query_end = data["min_start_query_end"]

        strand = determine_strand_for_pair(input_blast, data["max_end_query"],ref)
        original_max_end_query = data["max_end_query"]
        if strand == "-":
            data["max_end_query"],data["max_end_query_end"] =conver_minus_strand2_plus(data["max_end_query"], data["max_end_query_start"],fai_len)
        # print(strand, data["max_end_query"], data["max_end_query_end"])
        end_query = re.split('(\+|-)', data["max_end_query"])
        end_query = [end_query[n] + end_query[n + 1] for n in range(0, len(end_query) - 1, 2)]
        #end_query_start = data["max_end_query_start"]
        end_query_end = data["max_end_query_end"]

        start_query_filtered = []
        cumulative_len = 0
        for seg in start_query:
            seg_len = get_seg_len(seg, fai_len)
            current_pos = cumulative_len + seg_len
            fraction = float(current_pos - start_query_start)/float(seg_len)
            if cumulative_len + seg_len > start_query_start and fraction > 0.5:
                start_query_filtered.append(seg)
            cumulative_len += seg_len

        end_query_filtered = []
        cumulative_len = 0
        for seg in end_query:
            seg_len = get_seg_len(seg, fai_len)
            cumulative_len += seg_len
            fraction = float(cumulative_len - end_query_end)/float(seg_len)
            if cumulative_len < end_query_end or fraction < 0.5:
                end_query_filtered.append(seg)
        if data["min_start_query"] == data["max_end_query"]:
            intersection_list = [value for value in end_query_filtered if value in start_query_filtered]
            ref_start_end_segs[data["min_start_query"]]=intersection_list
            ref_start_end_segs[original_min_start_query] = intersection_list
        else:
            ref_start_end_segs[data["min_start_query"]]=start_query_filtered
            ref_start_end_segs[original_min_start_query] = start_query_filtered
            ref_start_end_segs[data["max_end_query"]] = end_query_filtered
            ref_start_end_segs[original_max_end_query] = end_query_filtered
        #print(f"Reference: {ref}")
        #print(f"Start Query: {start_query_filtered} {start_query_start}")
        #print(f"End Query: {end_query_filtered} {end_query_end}")
    return ref_start_end_segs
def get_seg_len(seg, fai_len):
    seg_p = seg.replace("+","").replace("-","").replace("\t","")
    return fai_len[seg_p]
def check_gene_or_score(line, genes,scores):
    vs = re.split(r'\+|-|\t',line)
    for v in vs:
        if v!= "":
            if v in genes.keys() or v in scores.keys():
                return True
    return False
def get_line_len(line, fai_len):
    result_len = 0
    vs = re.split(r'\+|-|\t',line)
    for v in vs:
        if v != "":
            result_len += get_seg_len(v, fai_len)
    #print(vs,line,"wewewe",result_len)
    return result_len
def main():
    parser = argparse.ArgumentParser(description="Filter by BLAST")
    parser.add_argument("input_file", help="Input file")
    parser.add_argument("cycle_txt", help="Cycle text")
    parser.add_argument("fasta_fai", help="Fasta fai")
    parser.add_argument("second_match", help="Second match file")
    parser.add_argument("run_model", help="Run model")
    parser.add_argument("blast_ratio", help="Blast ratio", type=float)
    parser.add_argument("blast_len_threshold", help="Blast length threshold", type=int)
    parser.add_argument("-s", "--single_ref", help="Single reference", default="")
    parser.add_argument("--gene_hit", help="Gene hit")
    parser.add_argument("--score", help="Score")
    parser.add_argument("--before_cut", help="contains info before cut(used for filter)")
    args = parser.parse_args()
    input_file = args.input_file
    cycle_txt = args.cycle_txt
    fasta_fai = args.fasta_fai
    second_match = open(args.second_match,"w")
    run_model = args.run_model
    blast_ratio = args.blast_ratio
    blast_len_threshold = args.blast_len_threshold
    single_ref = args.single_ref

    ref_list = {}
    genes={}
    scores={}
    with open(args.gene_hit, 'r') as gene_lst:
        for gene_r in gene_lst:
            contig_name,_=gene_r.strip().split("\t")
            genes[contig_name] = '1'

    with open(args.score, 'r') as score_lst:
        for score_r in score_lst:
            line_lst = score_r.strip().split("\t")
            scores[line_lst[0]] = str(line_lst[1])
    with open(args.input_file) as f:
        for line in f:
            line = line.strip("\n").split()
            if line[1] not in ref_list:
                ref_list[line[1]] = int(line[4])

    fai_len = {}
    with open(args.fasta_fai, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            fai_len[fields[0]] = int(fields[1])

    ref_contig = {}
    res = set()
    if args.run_model == "1":
        with open(args.cycle_txt) as r:
            for line in r.readlines():
                line_len = 0
                splited=re.split(r'[+-]',line.strip())
                for v in splited:
                    if v != "" or v !=" ":
                        line_len += get_line_len(v, fai_len)
                if line_len >= 10000:
                    liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","")
                    res.add(liner.strip("\n"))
# 生成需要第二步match的文件： format: res \t ref

    title_contig = {}
    ref_contig_l = {}
    blast_segs = set()
    prev_seg = ""
    prev_len = 0
    prev_ref = ""
    with open(args.input_file) as f_in:
        for line in f_in.readlines():

            t = line.strip().split("\t")
            #if "EDGE_28981966_length_1824_cov_10" in prev_seg:
            #    print(prev_len,"tdtd")
            if (args.single_ref != "" and  t[1] not in args.single_ref):
                continue
            if (prev_seg != t[0] and prev_seg != "") or (prev_ref != t[1] and prev_ref != ""):
            #    if "EDGE_28981966_length_1824_cov_10" in prev_seg:
            #        print(prev_len, elen,"xdxd")
                elen = get_line_len(prev_seg, fai_len)
                if float(prev_len) / float(elen) > blast_ratio or prev_len > blast_len_threshold or check_gene_or_score(t[0], genes,scores):
                #if float(prev_len) / float(elen) > blast_ratio:
                    #print(prev_seg, prev_len,elen,float(prev_len) / float(elen),blast_ratio,blast_len_threshold, "111111111")
                    blast_segs.add(prev_seg)
                prev_seg = t[0]
                prev_ref = t[1]
                prev_len = int(t[5])
            else:
                if float(t[2]) > 75:
            #       if "EDGE_28981966_length_1824_cov_10" in prev_seg:
            #           print(prev_len,float(t[2]),blast_ratio*100,"mmm")
                    prev_len = prev_len + int(t[5])
                prev_seg = t[0]
                prev_ref = t[1]
    elen = get_line_len(prev_seg, fai_len)
    if elen != 0:
        if float(prev_len) / float(elen) > blast_ratio or prev_len > blast_len_threshold:
        #if float(prev_len) / float(elen) > blast_ratio:
            #print(prev_seg, prev_len,elen,"22222222")
            blast_segs.add(t[0])
    ref_start_end_segs = cut_end_contig(args.input_file,blast_segs, fai_len,args.single_ref)
    for fline in open(args.input_file).readlines():
        # ref_length = ref_list[ref]
        contig_num = 10
        count = contig_num - 1
        contig = ""
        line = fline.strip("\n").split("\t")
        if args.single_ref != "" and line[1] not in args.single_ref:
            continue
        if line[0] not in blast_segs:
            continue
        #if (float(line[2]) < 75 or int(line[5])/int(line[3]) < 0.5) and int(line[5]) < 2000:
        #    continue
        if line[1] not in ref_contig.keys():
            ref_contig[line[1]] = []
            #title_contig[line[1]] = []
            ref_contig_l[line[1]] = 0
        if line[0] == contig:
            start = min(int(line[10]), int(line[11]))
            stop = max(int(line[10]), int(line[11]))
            # rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
            ref_contig[line[1]].append([start, stop, line[0]])
            #if line[0] not in title_contig[line[1]]:
            #    title_contig[line[1]].append(line[0])
            ref_contig_l[line[1]] = ref_contig_l[line[1]] + (stop - start)
            # currentAxis.add_patch(rect)
            # title.add(line[0])
        else:
            count -= 1
            contig = line[0]
            start = min(int(line[10]), int(line[11]))
            stop = max(int(line[10]), int(line[11]))
            # rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
            ref_contig[line[1]].append([start, stop, line[0]])
            #if line[0] not in title_contig[line[1]]:
            #    title_contig[line[1]].append(line[0])
            ref_contig_l[line[1]] = ref_contig_l[line[1]] + (stop - start)
            # currentAxis.add_patch(rect)
            # title.add(line[0])
    for key, value in ref_contig.items():
        title_contig[key] = []
        ref_contig[key] = sorted(value, key=itemgetter(1))
        for v in ref_contig[key]:
            if v[2] not in title_contig[key]:
                title_contig[key].append(v[2])
    contig_ref = {}
    for ref in ref_list:
        if ref not in ref_contig.keys():
            continue
        title = set()
        ref_length = ref_list[ref]
        cover = [0] * ref_length
        vs = ref_contig[ref]
        # if ref == "CP024485.1":
        for v in vs:
            # print("CP024485.1",v[0],v[1])
            # if v[-1] in contig_cov.keys():
            #     contig_cov[v[-1]] = contig_cov[v[-1]] + (v[1]-v[0])
            # else:
            #     contig_cov[v[-1]] = (v[1]-v[0])
            for i in range(v[0], v[1]):
                cover[i] = 1
        # print(cover)
        un_covered = cover.count(0)
        if ref in ref_contig.keys():
            # if (float(ref_contig_l[ref])/float(ref_length))<0.8:
            if un_covered / ref_length > 0.4:
                #print("not ok", ref, un_covered, un_covered / ref_length, ref_length)
                plt.close()
                continue
        pt = ""
        for i in title_contig[ref]:
            pt = pt +"\t"+ i
        if pt in contig_ref.keys():
            contig_ref[pt].append(ref)
        else:
            contig_ref[pt] = [ref]

    import re

    k_lens = {}
    for k in contig_ref.keys():
        k_lens[k] = []
        t_l = 0
        t = re.split(r"[+-]", k.strip())
        #print(t,"ewwe")
        for i in t:
            if i == "":
                continue
            l = get_line_len(i,fai_len)
            k_lens[k].append(l)

    # print("k_lens===============================================")
    # print(k_lens)
    # print("k_lens===============================================")

    import math

    count = 0
    result = []
    skip = []
    replace = {}
    similar_array = []
    for fk in k_lens.keys():
        if fk in skip:
            continue
        a = k_lens[fk]
        oflag = True
        for sk in k_lens.keys():
            b = k_lens[sk]
            if fk == sk or sk < fk or sk in skip:
                continue
            tmp = [j for j in a if j in b]
            if sum(tmp) / sum(a) > 0.8 or sum(tmp) / sum(b) > 0.8:
                oflag = False
                flag = True
                for suba in similar_array:
                    if fk in suba:
                        suba.append(sk)
                        flag = False
                        break
                    elif sk in suba:
                        suba.append(fk)
                        flag = False
                        break
                if flag:
                    similar_array.append([fk,sk])
                    # else:

                # min_diff_ref_with_contig = 1
                # min_ref = ""i
                # for ref in contig_ref[fk]:
                #     if math.abs(1-ref_list[ref]/sum(a)) >
                # pass
                # print(sum(tmp) / sum(a), sum(tmp) / sum(b))
                # count = count + 1
                # flag = False
                # print(fk,"----------",sk,"duplicate",tmp)
                # replace[fk] = sk
                # replace[sk] = fk
                # if sum(b) < sum(a):
                #     if sk not in result:
                #         result.append(sk)
                #     skip.append(fk)
                # else:
                #     if fk not in result:
                #         result.append(fk)
                #     skip.append(sk)
            # else:

        if oflag:
            similar_array.append([fk])

    # print(count)
    # # print(tmp)
    # print(len(contig_ref))
    for s in similar_array:
        max_v = 0
        max_it = ""
        for it in s:
            a = k_lens[it]
            if sum(a) > max_v:
                max_v = sum(a)
                max_it = it
        result.append(max_it)


    t = 0
    visited_path = []
    for k in result:
        if "self" in k or "gene" in k:
            pass
            # print(k)
        t = t + len(contig_ref[k])
        for ref in contig_ref[k]:
            # pass

            ref_length = ref_list[ref]
            contig_num = 30
            cover = list(range(0, ref_length))
            # plot
            plt.figure(figsize=(20, 10))
            plt.tight_layout()
            plt.axis('off')
            plt.xlim(xmin=0 - 5)
            plt.xlim(xmax=ref_length + 1000)
            plt.ylim(ymin=-0.5)
            plt.ylim(ymax=contig_num + 2)

            currentAxis = plt.gca()
            rect = patches.Rectangle((1, contig_num), width=ref_length, height=0.3, linewidth=2, edgecolor="black",
                                    facecolor="#549FCF")
            currentAxis.add_patch(rect)
            contig = ""
            count = contig_num - 1
            vs = ref_contig[ref]
            # if ref == "CP077390.1":
            # for v in vs:
            #     print(v[0],v[1])
            #     for i in range(v[0],v[1]):
            #         cover[i] = 1
            # un_covered = cover.count(0)
            # print(ref,ref_contig_l[ref],ref_length,float(ref_contig_l[ref])/float(ref_length))
            for v in vs:
                if (v[1] - v[0]) / ref_length < 0.05 and (v[1] - v[0]) < 300:
                    continue
                if v[-1] == contig:
                    rect = patches.Rectangle((v[0], count), width=(v[1] - v[0]), height=0.3, linewidth=2, edgecolor="black",
                                            facecolor="#FF8C00")
                    currentAxis.add_patch(rect)
                else:
                    count -= 1
                    contig = v[-1]
                    rect = patches.Rectangle((v[0], count), width=(v[1] - v[0]), height=0.3, linewidth=2, edgecolor="black",
                                            facecolor="#FF8C00")
                    currentAxis.add_patch(rect)
            # plt.title(pt)
            if k in replace.keys():
                k2 = replace[k]
            else:
                k2 = k
            if k2 not in visited_path:
                path = k2.replace("gene_score", "").replace("score", "").replace("gene", "").replace("self", "").replace("self-gene",                                                                                          "").replace("ref","")
                second_match.write(path.replace("\t","") + '\t' + ref + '\n')
                res.add(path.strip("\n"))
                plt.savefig(os.path.join(os.path.dirname(args.input_file), "N" + k2[0:15] + "_" + k2[-15:-1] + "%s_blast.png" % ref), dpi=300)
            visited_path.append(k2)
            plt.close()
    
    if args.before_cut:
        with open(args.before_cut,"w") as bc_fout:
            for item in res:
                new_item = ""
                for seg in item.strip().split("\t"):
                    if seg in ref_start_end_segs.keys():
                        seg = "".join(ref_start_end_segs[seg])
                    new_item += seg
                new_item_str=new_item.replace("\t","").replace("+","+\t").replace("-","-\t")
                print(new_item_str.strip())
                bc_fout.write(new_item_str+":"+item.replace("\t","").replace("+","+\t").replace("-","-\t")+"\n")
    second_match.close()

if __name__ == "__main__":
    main()
