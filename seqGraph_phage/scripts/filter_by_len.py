import sys

if len(sys.argv) != 2:
    print("Usage: python script.py <blast_output_file>")
    sys.exit(1)

input_file = sys.argv[1]

# Function to merge overlapping intervals
def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged_intervals = [sorted_intervals[0]]

    for current_interval in sorted_intervals[1:]:
        last_interval = merged_intervals[-1]

        if current_interval[0] <= last_interval[1]:
            new_interval = (last_interval[0], max(last_interval[1], current_interval[1]))
            merged_intervals[-1] = new_interval
        else:
            merged_intervals.append(current_interval)

    return merged_intervals

# Initialize a dictionary to store aligned subject references
aligned_subjects = {}

with open(input_file, "r") as file_A:
    for line in file_A:
        columns = line.strip().split("\t")

        query_id = columns[0]
        subject_id = columns[1]
        identity = float(columns[2])
        q_len = int(columns[3])
        s_len = int(columns[4])
        s_start = int(columns[8])
        s_end = int(columns[9])

        if identity > 95:
            if query_id not in aligned_subjects:
                aligned_subjects[query_id] = {subject_id: {'intervals': [(s_start, s_end)], 'q_len': q_len, 's_len': s_len}}
            else:
                if subject_id not in aligned_subjects[query_id]:
                    aligned_subjects[query_id][subject_id] = {'intervals': [(s_start, s_end)], 'q_len': q_len, 's_len': s_len}
                else:
                    aligned_subjects[query_id][subject_id]['intervals'].append((s_start, s_end))

# Merge overlapping intervals and calculate non-overlapping length
non_overlap_length = {}

for query_id, subjects in aligned_subjects.items():
    non_overlap_length[query_id] = {}

    for subject_id, data in subjects.items():
        merged_intervals = merge_intervals(data['intervals'])
        length = sum(abs(interval[1] - interval[0]) + 1 for interval in merged_intervals)
        non_overlap_length[query_id][subject_id] = {'length': length, 'q_len': data['q_len'], 's_len': data['s_len']}

# Print non-overlapping length for each query sequence and aligned subject reference
for query_id, subjects in non_overlap_length.items():
    for subject_id, data in subjects.items():
        print(f"Query: {query_id} (Length: {data['q_len']}), Subject: {subject_id} (Length: {data['s_len']}), Non-overlapping Length: {data['length']}")
