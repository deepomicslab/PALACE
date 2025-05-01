import re
import argparse
MIN_STRICT_GENE = 5
min_len_gene_dict = {2000:2,5000:5,10000:10,20000:20,50000:30,100000:35,200000:40,500000:50}
#min_len_gene_dict={2000:2,5000:5,10000:8,20000:12,100000:20,200000:35,500000:50}

#def check_gene(length, gene_count, min_len_gene_dict={2000:2,5000:5,10000:8,20000:12,100000:20,200000:35,500000:50}):
#    """
#    Check if the gene_count meets the minimum requirement for the given length.
#
#    :param length: The length of the sequence.
#    :param gene_count: The number of genes.
#    :param min_len_gene_dict: A dictionary mapping minimum length thresholds to required gene counts.
#    :return: True if gene_count meets or exceeds the required minimum, otherwise False.
#    """
#    # Sort the dictionary keys in ascending order
#    sorted_keys = sorted(min_len_gene_dict.keys())
#
#    # Find the highest applicable threshold
#    required_genes = 0
#    for key in sorted_keys:
#        if int(length) >= key:
#            required_genes = min_len_gene_dict[key]
#        else:
#            break
#
#    # Check if the gene count is sufficient
#    return gene_count >= required_genes
def check_gene(length, gene_count, min_gene_density=1):
    """
    Check if the gene_count meets the minimum required gene density for the given sequence length.

    :param length: The length of the sequence (bp).
    :param gene_count: The number of genes.
    :param min_gene_density: The minimum gene density (genes per 1000 bp).
    :return: True if gene_count meets or exceeds the required minimum, otherwise False.
    """
    if gene_count >= 40:
        return True
    else:
        # 计算最低需要的基因数
        required_genes = min_gene_density * (length / 1000)

        # 判断是否满足最低基因密度
    return gene_count >= required_genes

# Example usage:
#print(check_gene(2500, 3))  # Expected: True (2500 > 2000, requires >2 genes)
#print(check_gene(6000, 4))  # Expected: False (6000 > 5000, requires >5 genes)
#print(check_gene(15000, 9)) # Expected: True (15000 > 10000, requires >8 genes)
#print(check_gene(21000, 12))# Expected: True (21000 > 20000, requires >12 genes)
def parse_graph(file_path,gene_res):
    in_gene = []
    in_score = []
    both = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.split()
            if columns[0] == 'SEG':
                try:
                    fourth_value = float(columns[4])
                    fifth_value = float(columns[5])
                    if columns[1] in gene_res and fifth_value > 0.7:
                        both.append(columns[1])
                    elif fourth_value > 0.9:
                        in_gene.append(columns[1])
                    elif fifth_value > 0.7:
                        in_score.append(columns[1])
                except (IndexError, ValueError):
                    # Handle lines that don't have enough columns or have non-float values
                    continue
    return in_gene,in_score,both

def get_edge_len(edge):
    return int(edge.split("_")[3])
def get_path_len(paths):
    full_len = 0
    for p in paths:
        if len(p) == 0 or p == "+" or p == "-" or p == " ":
            continue
        pl = int(p.split("_")[3])
        full_len += pl
    return full_len
def parse_result(file_path):
    result = []
    pattern = re.compile(r'\t+')
    with open(file_path, 'r') as file:
        for line in file:
            if len(line.strip()) == 0:
                continue
            if "iter" in line:
                continue
            line = line.replace("+","+\t").replace("-","-\t")
            items = pattern.split(line.strip())
            result.append(items)
    return result


def get_items_in_keeped(items, in_gene,in_score,in_both,strict=dict()):
    gene_score_dict = []
    total_gene_count = 0
    gene_len = 0
    score_len = 0
    both_len = 0
    for tmp_item in items:
        item = tmp_item.replace("+","").replace("-","").replace(" ","").replace("\t","")
        if item in strict.keys():
            total_gene_count += int(strict[item])
        if item in in_both:
            gene_score_dict.append((tmp_item,2))
            both_len += get_edge_len(item)
        elif item in strict.keys():
            if len(strict) != 0:
                if check_gene(get_edge_len(item),strict[item]):
                    gene_score_dict.append((tmp_item, 1))
                    gene_len += get_edge_len(item)
                else:
                    gene_score_dict.append((tmp_item, -1))
            else:
                gene_score_dict.append((tmp_item, 1))
                gene_len += get_edge_len(item)
        elif item in in_score:
            gene_score_dict.append((tmp_item, 0))
            score_len += get_edge_len(item)
        else:
            gene_score_dict.append((tmp_item, -1))


    return float(gene_len),float(score_len),float(both_len),gene_score_dict,total_gene_count

def split_list(arr):
    """
    Splits `arr` into sublists according to the rule:
    1. Traverse the list from left to right.
    2. Accumulate items in a current sublist.
    3. Whenever you encounter consecutive items with value == -1:
       - Compute the sum of get_len(key) for those items.
       - If the sum >= 2000, remove those -1 items entirely,
         finalize (close) the current sublist right before them,
         and start a new sublist after them.
       - Otherwise, include them in the current sublist.
    4. Finally, flatten all sublists into a single list containing only the keys.
    5. For any blocks of -1 that caused a split (sum >= 2000),
       those -1 items are not included at all in the final result.
    """
    result_sublists = []
    current_sublist = []
    i = 0
    n = len(arr)

    while i < n:
        key, val = arr[i]

        if val != -1:
            # Just append normal items (value != -1)
            current_sublist.append((key, val))
            i += 1
        else:
            # Found a block of consecutive -1 items
            j = i
            block_len_sum = 0
            while j < n and arr[j][1] == -1:
                block_len_sum += get_edge_len(arr[j][0])
                j += 1

            # block_len_sum is the total length for this consecutive -1 block
            if block_len_sum >= 1000:
                # Remove the entire block: do not include in current_sublist
                # Finalize and close the current_sublist right before this block
                if current_sublist:
                    result_sublists.append(current_sublist)
                current_sublist = []
            else:
                # If below the threshold, we keep these -1 items
                while i < j:
                    current_sublist.append(arr[i])
                    i += 1

            # Move on past the block of -1s if not already moved
            i = j

    # Finalize whatever remains
    if current_sublist:
        result_sublists.append(current_sublist)

    # Now flatten while keeping only the keys
    # Note: any block of -1s that caused a split has already been discarded,
    # so they never appear in the final result anyway.
    final_result = []
    for sublist in result_sublists:
        keys = []
        for (key,val) in sublist:
            keys.append(key)
        final_result.append(keys)
    return final_result

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process remain results")

    # Adding arguments
    parser.add_argument("graph", help="The file path to the graph file.")
    parser.add_argument("remain", help="The file for processing.")
    parser.add_argument("output", help="output.")
    parser.add_argument("threshold", type=float, help="the length percent to keeped")
    parser.add_argument("minlen", type=float, help="the length percent to keeped")
    parser.add_argument("beforcut", help="cutted file, this file not effect result")
    parser.add_argument("gene_file", help="gene_file")

    # Parse arguments
    args = parser.parse_args()
    gene_res = dict()
    with open(args.gene_file, 'r') as gene_lst:
        for gene_r in gene_lst:
            contig_name, gene_count = gene_r.split("\t")
            gene_res[contig_name] = int(gene_count)

    in_gene,in_score,in_both = parse_graph(args.graph,gene_res)
    results1 = parse_result(args.remain)
    final_keeped = []
    #TODO, this invoke get_items_in_keeped 3 times, refine
    for items in results1:
        gene_len,score_len,both_len,gene_score_dict,total_gene = get_items_in_keeped(items, in_gene,in_score,in_both, gene_res)
        len2 = float(get_path_len(items))
        if len2 < args.minlen:
            continue
        print(items,gene_len,score_len,both_len,len2,both_len/len2 >= args.threshold/2,(gene_len + score_len + both_len)/len2 >= args.threshold,len2)
        if (both_len/len2 >= args.threshold/2 and (gene_len + score_len + both_len)/len2 >= args.threshold) or gene_len == len2 and len2 >= args.minlen:
            final_keeped.append(items)
        else:
            _,_,_,gene_score_dict,_ = get_items_in_keeped(items, in_gene,in_score,in_both,gene_res)
            sublists = split_list(gene_score_dict)
            for idx, sublst in enumerate(sublists, start=1):
                gene_len,score_len,both_len,gene_score_dict,total_gene = get_items_in_keeped(sublst, in_gene,in_score,in_both, gene_res)
                len2 = float(get_path_len(sublst))
                print(total_gene,sublst,gene_len,score_len,both_len,len2,both_len/len2 >= args.threshold/2,(gene_len + score_len + both_len)/len2 >= args.threshold,len2)
                len2 = float(get_path_len(sublst))
                if ((float(gene_len)/len2 > 0.95) or float(gene_len + both_len)/len2 > 0.95 or float(both_len)/len2 > 0.95) and len2 >= args.minlen and total_gene >=8:
                    final_keeped.append(sublst)
    with open(args.output, "w") as f:
        for items in final_keeped:
            f.write("\t".join(items)+"\n")
    with open (args.beforcut, "w") as f:
        for items in final_keeped:
            f.write("\t".join(items)+":"+"\t".join(items)+"\n")
