import os
import sys
from collections import defaultdict
import re

def reverse_sign(char):
    """
    Reverses the sign of a character.
    This is a placeholder function. Replace it with the actual logic for reversing the sign.
    """
    if char == '+':
        return '-'
    elif char == '-':
        return '+'
    else:
        return char

def reverse_string(s):
    # Split the string by '+' or '-' and keep the delimiters
    parts = re.split(r'(\+|-)', s)

    # Combine each part with its following delimiter
    combined_parts = [parts[i] + parts[i+1] for i in range(0, len(parts) - 1, 2)]
    
    # Reverse the order of the combined parts
    combined_parts.reverse()

    # Reverse the sign of the last character in each part
    for i in range(len(combined_parts)):
        if combined_parts[i]:
            combined_parts[i] = combined_parts[i][:-1] + reverse_sign(combined_parts[i][-1])

    # Join the parts back together
    reversed_str = ''.join(combined_parts)

    return reversed_str
def read_file_content(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def process_line(directory, line):
    # Dictionary to count how many times each content appears
    content_count = defaultdict(int)
    
    # Split the line by comma to get the reference names
    refs = line.strip().split(',')
    
    for ref in refs:
        ref = ref.replace("|","_")
        ragtag_file = os.path.join(directory, f"{ref}_ragtag_scaffold_part.txt")
        
        # Read the content of the corresponding file
        if os.path.isfile(ragtag_file):
            content = read_file_content(ragtag_file)
            content.strip()
            if content in content_count.keys():
                content_count[content] += 1
            elif reverse_string(content) in content_count.keys():
                content_count[reverse_string(content)] += 1
            else:
                content_count[content] = 1

        else:
            print(f"Warning: File {ragtag_file} not found.")
    
    # Find the content that appears most frequently
    if content_count:
        most_frequent_content = max(content_count, key=content_count.get)
        return most_frequent_content
    else:
        return None

def main(directory, input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'a') as outfile:
        for line in infile:
            most_frequent_content = process_line(directory, line)
            if most_frequent_content:
                outfile.write(most_frequent_content+"\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <directory> <input_file> <output_file>")
        sys.exit(1)
    
    dir_path = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(dir_path, input_file, output_file)
