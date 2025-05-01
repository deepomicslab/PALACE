import argparse

def remove_duplicate_pairs(input_file, output_file):
    # Read the input file
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Check if the number of lines is even, if not, add an empty line
    if len(lines) % 2 != 0:
        lines.append('\n')
    
    # Group lines into pairs
    pairs = [(lines[i], lines[i + 1]) for i in range(0, len(lines), 2)]
    
    # Remove duplicate pairs
    unique_pairs = []
    seen_pairs = set()
    
    for pair in pairs:
        if pair not in seen_pairs:
            unique_pairs.append(pair)
            seen_pairs.add(pair)
    
    # Write the unique pairs to the output file
    with open(output_file, 'w') as file:
        for pair in unique_pairs:
            file.write(pair[0])
            file.write(pair[1])
    
    print(f"Unique pairs have been written to {output_file}")

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Remove duplicate pairs of lines from a file.")
    
    # Add arguments
    parser.add_argument('input_file', help="Path to the input file")
    parser.add_argument('output_file', help="Path to the output file")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function with command line arguments
    remove_duplicate_pairs(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
