import sys
import re

def extract_edges_from_file_b(filename):
    """Extract all EDGE identifiers from file b (ignoring +/- signs)"""
    edges = set()
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read()
            # Use regex to match EDGE pattern
            edge_pattern = r'EDGE_\d+_length_\d+_cov_[\d.]+'
            matches = re.findall(edge_pattern, content)
            for match in matches:
                edges.add(match)
    except FileNotFoundError:
        print(f"Error: File not found {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)
    
    return edges

def filter_lines_from_file_a(filename_a, edges_to_remove):
    """Filter lines in file a, remove lines containing any specified EDGEs"""
    filtered_lines = []
    
    try:
        with open(filename_a, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Check if this line contains any EDGE to be removed
                should_remove = False
                edge_pattern = r'EDGE_\d+_length_\d+_cov_[\d.]+'
                matches = re.findall(edge_pattern, line)
                
                for match in matches:
                    if match in edges_to_remove:
                        should_remove = True
                        break
                
                if not should_remove:
                    filtered_lines.append(line)
                    
    except FileNotFoundError:
        print(f"Error: File not found {filename_a}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {filename_a}: {e}")
        sys.exit(1)
    
    return filtered_lines

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <file_a> <file_b> <output_file>")
        print("Example: python script.py file_a.txt file_b.txt output.txt")
        sys.exit(1)
    
    file_a = sys.argv[1]
    file_b = sys.argv[2]
    output_file = sys.argv[3]
    
    print(f"Processing files: {file_a} and {file_b}")
    
    # Extract all EDGE identifiers from file b (ignoring +/- signs)
    print("Extracting EDGE identifiers from file b...")
    edges_in_b = extract_edges_from_file_b(file_b)
    print(f"Found {len(edges_in_b)} unique EDGEs in file b")
    
    # Filter lines in file a
    print("Filtering file a...")
    filtered_lines = filter_lines_from_file_a(file_a, edges_in_b)
    
    # Write results to specified output file
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            for line in filtered_lines:
                f.write(line + '\n')
        
        print(f"Filtering complete! Results saved to: {output_file}")
        print(f"Retained {len(filtered_lines)} lines")
        
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
