def reverse_and_flip_directions(input_string):
    """
    Reverse the sequence of segments in the input string while flipping their orientations.
    Each segment is presumed to end with a '+' or '-' indicating its orientation.

    :param input_string: str, the input string with segments ending in '+' or '-'
    :return: str, the reversed string with flipped orientations
    """
    # Split the string into segments based on the end characters '+' and '-'
    segments = []
    current_segment = ""

    for char in input_string:
        if char in '+-':
            # Add the current segment with the orientation character
            segments.append(current_segment + char)
            current_segment = ""
        else:
            current_segment += char

    # Reverse the list of segments
    segments.reverse()

    # Flip the orientation of each segment
    flipped_segments = []
    for segment in segments:
        if segment[-1] == '+':
            flipped_segments.append(segment[:-1] + '-')
        else:
            flipped_segments.append(segment[:-1] + '+')

    # Join the segments back into a single string
    return ''.join(flipped_segments)
def process_contig_line(contig_line):
    """
    Process the contig line to reverse and flip orientations based on the rules given.
    """
    contigs = contig_line.split('+')
    
    # Remove empty strings if any (due to trailing +)
    contigs = [contig for contig in contigs if contig.strip()]
    
    # Reverse the list of contigs and flip the orientation
    processed_contigs = []
    for contig in contigs:
        print(contig)
        if contig.endswith('-'):
            new_contig = contig[:-1] + '+'
        elif contig.endswith('+'):
            new_contig = contig[:-1] + '-'
        else:
            new_contig = contig + '-'  # No orientation change if no + or -
        processed_contigs.append(new_contig)
    
    # Join back with '+' and return
    return '+'.join(reversed(processed_contigs))

def process_file(input_filename, output_filename, is_remain):
    """
    Process the file to modify specific lines based on given conditions and write the output to a new file.
    """
    if is_remain:
        preref = ""
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                # Check if the line meets the specific conditions
                if len(columns) >= 9 and columns[0].endswith('_RagTag') and columns[4] == 'W':
                    if preref != columns[0] and preref != "":
                        print(preref)
                        outfile.write("\n")
                    if columns[8] == '-':
                        columns[5] = reverse_and_flip_directions(columns[5])
                        # Process the contig line to reverse and flip
                    
                    # Write the modified sixth column to the output file
                    outfile.write(columns[5])
                    preref =columns[0]
                elif columns[4] == 'W':
                    outfile.write(columns[5])
                    outfile.write("\n")
    else:
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line in infile:
                columns = line.strip().split()
                # Check if the line meets the specific conditions
                if len(columns) >= 9 and columns[0].endswith('_RagTag') and columns[4] == 'W':
                    if columns[8] == '-':
                        columns[5] = reverse_and_flip_directions(columns[5])
                        # Process the contig line to reverse and flip
                    
                    # Write the modified sixth column to the output file
                    outfile.write(columns[5])
            outfile.write("\n")

# Example usage
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_filename> <output_filename> <is_remain 0 or 1>")
    else:
        process_file(sys.argv[1], sys.argv[2], sys.argv[3] == "1")
