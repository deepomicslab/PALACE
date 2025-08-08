import argparse
from collections import defaultdict

def parse_blast_file(filename, ref_queries_output):
    ref_queries = defaultdict(list)
    query_ref_lengths = defaultdict(lambda: defaultdict(int))
    query_lengths = {}

    # 读取并解析 BLAST 输出文件
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 14:
                continue
            query_id = parts[0]      # qaccver
            ref_id = parts[1]        # saccver
            query_length = int(parts[3])  # qlen
            aligned_length = int(parts[5])  # length
            # 避免比对太琐碎了
            if aligned_length < 100 and float(aligned_length)/float(query_length) < 0.05:
                continue

            # 存储查询长度
            query_lengths[query_id] = query_length

            # 记录每个查询-参考对的比对长度
            query_ref_lengths[query_id][ref_id] += aligned_length

    # 根据条件过滤记录
    for query_id, ref_lengths in query_ref_lengths.items():
        for ref_id, total_aligned_length in ref_lengths.items():
            query_length = query_lengths[query_id]
            #if total_aligned_length >= 2000 or (total_aligned_length / query_length) >= 0.7:
            if (total_aligned_length / query_length) >= 0.7:
                ref_queries[ref_id].append(query_id)

    # 查找重叠的参考并合并成独立的集合
    alignments = [(query_id, ref_id, aligned_length) for query_id, ref_lengths in query_ref_lengths.items() for ref_id, aligned_length in ref_lengths.items()]
    num_alignments = len(alignments)

    ref_sets = []
    ref_to_set = {}

    for i in range(num_alignments):
        for j in range(i + 1, num_alignments):
            query_a, ref_a, len_a = alignments[i]
            query_b, ref_b, len_b = alignments[j]

            if ref_a != ref_b and query_a == query_b:
                overlap = min(len_a, len_b)
                if overlap / len_a > 0.7 or overlap / len_b > 0.7:
                    if ref_a in ref_to_set and ref_b in ref_to_set:
                        if ref_to_set[ref_a] != ref_to_set[ref_b]:
                            set_a = ref_to_set[ref_a]
                            set_b = ref_to_set[ref_b]
                            new_set = set_a.union(set_b)
                            for ref in new_set:
                                ref_to_set[ref] = new_set
                            ref_sets.remove(set_a)
                            ref_sets.remove(set_b)
                            ref_sets.append(new_set)
                    elif ref_a in ref_to_set:
                        ref_to_set[ref_b] = ref_to_set[ref_a]
                        ref_to_set[ref_a].add(ref_b)
                    elif ref_b in ref_to_set:
                        ref_to_set[ref_a] = ref_to_set[ref_b]
                        ref_to_set[ref_b].add(ref_a)
                    else:
                        new_set = set([ref_a, ref_b])
                        ref_to_set[ref_a] = new_set
                        ref_to_set[ref_b] = new_set
                        ref_sets.append(new_set)

    # 输出 ref_queries 到文件
    with open(ref_queries_output, 'w') as file:
        for ref, queries in ref_queries.items():
            file.write(f"{''.join(queries)}\t{ref}\n")

# 解析命令行参数
def main():
    parser = argparse.ArgumentParser(description='Parse BLAST output and find overlapping references.')
    parser.add_argument('input_file', help='Path to the BLAST output file')
    parser.add_argument('ref_queries_output', help='Path to the output file for reference queries')

    args = parser.parse_args()

    parse_blast_file(args.input_file, args.ref_queries_output)

if __name__ == "__main__":
    main()
