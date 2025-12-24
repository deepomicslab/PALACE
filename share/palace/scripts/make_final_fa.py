"""
Generate final FASTA sequences from processed sequence paths.

This script reads sequence path information and assembles final FASTA sequences
with proper orientation handling and custom prefixes.
"""

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq


class FinalFastaGenerator:
    """Generate final FASTA sequences from sequence paths."""
    
    def __init__(self, fasta_file, paths_file, output_file, sequence_prefix):
        """
        Initialize the FASTA generator.
        
        Args:
            fasta_file (str): Input FASTA file with source sequences
            paths_file (str): File containing sequence paths
            output_file (str): Output FASTA file path
            sequence_prefix (str): Prefix for output sequence names
        """
        self.fasta_file = fasta_file
        self.paths_file = paths_file
        self.output_file = output_file
        self.sequence_prefix = sequence_prefix
        
        self.sequence_records = {}
        self.sequence_counter = 0
    
    def load_sequence_records(self):
        """Load sequence records from the input FASTA file."""
        self.sequence_records = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))
    
    def clean_segment_name(self, segment):
        """
        Clean segment name by removing prefixes.
        
        Args:
            segment (str): Original segment name
            
        Returns:
            str: Cleaned segment name
        """
        return segment.replace("ref", "")
    
    def get_sequence_with_orientation(self, segment):
        """
        Get sequence with proper orientation handling.
        
        Args:
            segment (str): Segment name with orientation marker
            
        Returns:
            Bio.Seq.Seq: Sequence with correct orientation
        """
        if not segment:
            return Seq("")
        
        # Clean the segment name
        cleaned_segment = self.clean_segment_name(segment)
        
        # Determine sequence name and orientation
        if cleaned_segment.endswith(('+', '-')):
            sequence_name = cleaned_segment[:-1]
            orientation = cleaned_segment[-1]
        else:
            sequence_name = cleaned_segment
            orientation = '+'
        
        # Get the base sequence
        if sequence_name not in self.sequence_records:
            print(f"Warning: Sequence '{sequence_name}' not found in FASTA file")
            return Seq("")
        
        base_sequence = self.sequence_records[sequence_name].seq
        
        # Apply orientation
        if orientation == '-':
            return base_sequence.reverse_complement()
        else:
            return base_sequence
    
    def assemble_sequence_from_path(self, path_segments):
        """
        Assemble a complete sequence from path segments.
        
        Args:
            path_segments (list): List of segment names with orientations
            
        Returns:
            Bio.Seq.Seq: Assembled sequence
        """
        assembled_sequence = Seq("")
        
        for segment in path_segments:
            if segment:  # Skip empty segments
                segment_sequence = self.get_sequence_with_orientation(segment)
                assembled_sequence += segment_sequence
        
        return assembled_sequence
    
    def generate_sequence_name(self):
        """
        Generate a unique sequence name for output.
        
        Returns:
            str: Generated sequence name
        """
        self.sequence_counter += 1
        return f"{self.sequence_prefix}_phage_{self.sequence_counter}"
    
    def should_skip_line(self, line):
        """
        Check if a line should be skipped during processing.
        
        Args:
            line (str): Input line to check
            
        Returns:
            bool: True if line should be skipped
        """
        skip_keywords = ["all", "iter", "loop"]
        return any(keyword in line for keyword in skip_keywords)
    
    def process_paths_and_generate_fasta(self):
        """Process path file and generate final FASTA sequences."""
        with open(self.paths_file, 'r') as paths_input, \
             open(self.output_file, 'w') as fasta_output:
            
            for line in paths_input:
                line = line.strip()
                
                # Skip header and control lines
                if self.should_skip_line(line):
                    continue
                
                # Parse segments from the line
                path_segments = re.split(r'\t+', line.strip('\n'))
                path_segments = [seg for seg in path_segments if seg]  # Remove empty segments
                
                if not path_segments:
                    continue
                
                # Assemble sequence from path
                assembled_sequence = self.assemble_sequence_from_path(path_segments)
                
                if len(assembled_sequence) > 0:
                    # Generate output sequence name
                    sequence_name = self.generate_sequence_name()
                    
                    # Write to output FASTA
                    fasta_output.write(f">{sequence_name}\n")
                    fasta_output.write(f"{str(assembled_sequence)}\n")
    
    def run(self):
        """Execute the complete FASTA generation pipeline."""
        print("Loading sequence records...")
        self.load_sequence_records()
        
        print("Processing paths and generating FASTA...")
        self.process_paths_and_generate_fasta()
        
        print(f"Generated {self.sequence_counter} sequences in {self.output_file}")


def main():
    """Main function to handle command line arguments and run generation."""
    if len(sys.argv) != 5:
        print("Usage: python make_final_fa.py <fasta_file> <paths_file> "
              "<output_file> <sequence_prefix>")
        sys.exit(1)
    
    generator = FinalFastaGenerator(
        fasta_file=sys.argv[1],
        paths_file=sys.argv[2],
        output_file=sys.argv[3],
        sequence_prefix=sys.argv[4]
    )
    
    generator.run()


if __name__ == "__main__":
    main()