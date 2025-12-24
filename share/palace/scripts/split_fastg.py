#!/usr/bin/env python
##This scripts copied from https://github.com/Shamir-Lab/Recycler/blob/master/recyclelib/utils.py
import re, argparse, os

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'make_fasta_from_fastg converts fastg assembly graph to fasta format'
    )
    parser.add_argument('-g','--graph',
                        help='(spades 3.50+) FASTG file to process [recommended: before_rr.fastg]',
                        required=True, type=str
                        )
    parser.add_argument('-o','--output',
                        help='output file name for FASTA of cycles',
                        required=False, type=str
                        )

    return parser.parse_args()

def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader: 
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def parse_lines(fastg, ofile):
    fp = open(fastg, 'r')
    seen = set() ##
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1]).split(" ")[0]
        if name[-1] == "'":
            name = name[:-1]
            seq = reverse_complement(seq)
        if name in seen: continue
        else: seen.add(name)
        line = ">"+name+"\n"+seq+"\n"
        ofile.write(line)
def reverse_complement(dna):
    """
    This function takes a DNA sequence as input, converts it to uppercase,
    and returns its reverse complement.

    Parameters:
    dna (str): A string representing the DNA sequence.

    Returns:
    str: The reverse complement of the given DNA sequence in uppercase.
    """
    # Convert the input DNA sequence to uppercase
    dna = dna.upper()

    # Define the complement mapping for each nucleotide
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    # Generate the complement sequence
    complement_sequence = ''.join(complement[base] for base in dna)

    # Reverse the complement sequence to get the reverse complement
    reverse_complement_sequence = complement_sequence[::-1]

    return reverse_complement_sequence
if __name__=='__main__':
    args = parse_user_input()
    fastg = args.graph
    fp = open(fastg, 'r')
    files_dir = os.path.dirname(fp.name)

    # output 1 - fasta of sequences
    if args.output:
        fasta_ofile = args.output
    else:
        (root,ext) = os.path.splitext(fp.name)
        fasta_ofile = root + ext.replace(".fastg", ".nodes.fasta")

    f_nodes_fasta = open(fasta_ofile, 'w')
    parse_lines(fastg, f_nodes_fasta)
