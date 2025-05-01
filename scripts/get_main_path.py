import re
import sys

def parse_file_fa(file_path):
    relevant_items = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("SEG"):
                parts = line.split()
                if float(parts[-1]) > -2:
                    relevant_items.append(parts[1])
    return relevant_items

def get_path_len(path):
    sum = 0
    for item in path:
        if item.startswith("EDGE"):
            arr = item.split("_")
            sum +=int(arr[3])
    return sum
def analyze_file_fb(file_path, relevant_items):
    max_count = 0
    most_frequent_line = None
    result = []
    with open(file_path, 'r') as file:
        for line in file:
            items = re.split(r'\t+', line.strip())
            total_len = get_path_len(items)
            in_blast_items = [item for item in items if item[:-1] in relevant_items]
            count = len(in_blast_items)
            in_blast_len = get_path_len(in_blast_items)
            if float(in_blast_len)/float(total_len) >= 0.9 and in_blast_len > 2000:
                result.append(line.strip())
            if count > max_count:
                max_count = count
                most_frequent_line = line.strip()
    result.append(most_frequent_line)
    return result

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <file_fa_path> <file_fb_path> <output_file_path>")
        sys.exit(1)

    file_fa = sys.argv[1]
    file_fb = sys.argv[2]
    output_file_path = sys.argv[3]

    relevant_items_from_fa = parse_file_fa(file_fa)
    lines = analyze_file_fb(file_fb, relevant_items_from_fa)

    with open(output_file_path, 'w') as output_file:
        for line in lines:
            output_file.write(line + '\n')

    print(f"The line that appears most frequently has been written to '{output_file_path}'.")

if __name__ == "__main__":
    main()
