"""
FASTA sequence assembler from path descriptions.

This script reads contig paths from a file and assembles sequences from a FASTA file
based on the path information, handling forward and reverse complement orientations.
"""

import pysam
import sys


def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence.
    Handles both uppercase and lowercase nucleotides.
    
    Args:
        sequence (str): DNA sequence containing A/a, T/t, G/g, C/c
        
    Returns:
        str: Reverse complement of the input sequence, preserving case
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    return ''.join(complement_map.get(base, base) for base in reversed(sequence))


def fetch_contig_sequence(fasta_file, contig_name):
    """
    Fetch sequence for a contig, with fallback name handling.
    
    Args:
        fasta_file: pysam.FastaFile object
        contig_name (str): Name of the contig to fetch
        
    Returns:
        str: DNA sequence of the contig
    """
    try:
        return fasta_file.fetch(contig_name)
    except Exception:
        # Fallback: try removing the last underscore-separated part
        fallback_name = "_".join(contig_name.split('_')[:-1])
        return fasta_file.fetch(fallback_name)


def process_contig(fasta_file, contig_token):
    """
    Process a single contig token and return its sequence.
    
    Args:
        fasta_file: pysam.FastaFile object
        contig_token (str): Contig name with optional orientation (+ or -)
        
    Returns:
        str: Processed DNA sequence
    """
    contig_token = contig_token.replace(" ", "").strip()
    
    if len(contig_token) <= 1:
        return ""
    
    # Handle oriented contigs (ending with + or -)
    if contig_token.endswith('+'):
        contig_name = contig_token[:-1]
        if len(contig_name) == 0:
            print(f"Warning: Empty contig name for token '{contig_token}'")
            return ""
        return fetch_contig_sequence(fasta_file, contig_name)
    
    elif contig_token.endswith('-'):
        contig_name = contig_token[:-1]
        if len(contig_name) == 0:
            print(f"Warning: Empty contig name for token '{contig_token}'")
            return ""
        sequence = fetch_contig_sequence(fasta_file, contig_name)
        return reverse_complement(sequence)
    
    else:
        # No orientation specified, use forward strand
        try:
            return fetch_contig_sequence(fasta_file, contig_token)
        except Exception:
            fallback_name = "_".join(contig_token.split('_')[:-1])
            print(f"Contig not found: {fallback_name}")
            return fetch_contig_sequence(fasta_file, fallback_name)


def should_skip_line(line):
    """
    Check if a line should be skipped during processing.
    
    Args:
        line (str): Input line to check
        
    Returns:
        bool: True if line should be skipped
    """
    return (line.startswith("iter") or 
            line.startswith("self") or 
            line.strip() == "")


def write_fasta_record(output_file, header, sequence):
    """
    Write a FASTA record to the output file.
    
    Args:
        output_file: File handle for output
        header (str): FASTA header (without >)
        sequence (str): DNA sequence
    """
    output_file.write(f">{header}\n")
    output_file.write(f"{sequence}\n")


def main():
    """Main function to process path file and generate FASTA output."""
    if len(sys.argv) != 5:
        print("Usage: python make_fa_from_path.py <fasta_file> <paths_file> <output_file> <mode>")
        sys.exit(1)
    
    fasta_filename = sys.argv[1]
    paths_filename = sys.argv[2]
    output_filename = sys.argv[3]
    output_mode = sys.argv[4]
    
    print("make_fa_from_path.py running")
    
    # Open input files
    fasta_file = pysam.FastaFile(fasta_filename)
    
    try:
        with open(paths_filename, 'r') as paths_file, \
             open(output_filename, 'w') as output_file:
            
            for line_index, line in enumerate(paths_file):
                if should_skip_line(line):
                    continue
                
                # Parse contigs from the line
                contig_tokens = line.strip().split('\t')
                
                # Assemble sequence from all contigs in the path
                assembled_sequence = ""
                for contig_token in contig_tokens:
                    sequence_part = process_contig(fasta_file, contig_token)
                    assembled_sequence += sequence_part
                
                # Generate output header based on mode
                if output_mode == "0":
                    header = f"res_{line_index + 1}_{len(assembled_sequence)}"
                else:
                    header = "".join(contig_tokens)
                
                # Write the assembled sequence
                write_fasta_record(output_file, header, assembled_sequence)
    
    finally:
        fasta_file.close()


if __name__ == "__main__":
    main()