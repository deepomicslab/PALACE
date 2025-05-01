import sys
import re
from Bio import SeqIO
fa = open(sys.argv[1])
fai = open(sys.argv[2])
ref_file=sys.argv[3]
out_ref_fasta = open(sys.argv[4],"w")
out_ref_percent = open(sys.argv[5],"w")

def process_file(input_file):
    result = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("ref_index"):
                parts = line.split()
                first_int = None
                last_float = None

                # Find the first integer value after ref_index
                for part in parts:
                    if part.isdigit():
                        first_int = int(part)
                        break

                # Find the last float value in the line
                for part in reversed(parts):
                    try:
                        last_float = float(part)
                        break
                    except ValueError:
                        continue
            result[first_int] = last_float
    return result

fai_array = []
for line in fai.readlines():
    t = line.split("\t")
    fai_array.append(t[0])

record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

ref_index = process_file(ref_file)
for k in ref_index.keys():
    # t = re.split(r"\s+",line)
    ref_name = fai_array[int(k)-1]
    out_ref_fasta.write(">"+ref_name+"\n")
    out_ref_fasta.write(str(record_dict[ref_name].seq)+"\n")
    out_ref_percent.write(ref_name+"\t"+str(ref_index[k]) +"\n")