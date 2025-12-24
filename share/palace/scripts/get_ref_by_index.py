import sys
import re
from Bio import SeqIO


def parse_ref_file(ref_file):
    """Parse reference file to extract index and percentage values"""
    ref_data = {}
    
    with open(ref_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith("ref_index"):
                continue
            
            parts = line.split()
            
            # Find first integer (index)
            index = None
            for part in parts[1:]:  # Skip "ref_index"
                if part.isdigit():
                    index = int(part)
                    break
            
            # Find last float (percentage)
            percentage = None
            for part in reversed(parts):
                try:
                    percentage = float(part)
                    break
                except ValueError:
                    continue
            
            if index is not None and percentage is not None:
                ref_data[index] = percentage
    
    return ref_data


def load_fai_index(fai_file):
    """Load FAI index file into a dictionary"""
    index_to_name = {}
    
    with open(fai_file, 'r') as f:
        for idx, line in enumerate(f, 1):
            seq_name = line.split('\t')[0]
            index_to_name[idx] = seq_name
    
    return index_to_name


def main():
    """Main processing function"""
    # Parse command line arguments
    fasta_file = sys.argv[1]
    fai_file = sys.argv[2]
    ref_file = sys.argv[3]
    out_fasta = sys.argv[4]
    out_percent = sys.argv[5]
    
    # Load data
    print("Loading FASTA sequences...")
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    print("Loading FAI index...")
    index_to_name = load_fai_index(fai_file)
    
    print("Processing reference file...")
    ref_data = parse_ref_file(ref_file)
    
    # Write output files
    print("Writing output files...")
    with open(out_fasta, 'w') as fasta_out, open(out_percent, 'w') as percent_out:
        for index, percentage in sorted(ref_data.items()):
            # Get sequence name from index
            if index in index_to_name:
                seq_name = index_to_name[index]
                
                # Write FASTA sequence
                if seq_name in record_dict:
                    fasta_out.write(f">{seq_name}\n")
                    fasta_out.write(f"{record_dict[seq_name].seq}\n")
                    
                    # Write percentage
                    percent_out.write(f"{seq_name}\t{percentage}\n")
                else:
                    print(f"Warning: Sequence '{seq_name}' not found in FASTA file")
            else:
                print(f"Warning: Index {index} not found in FAI file")
    
    print("Processing complete!")


if __name__ == "__main__":
    main()
