import pysam
import sys


def rev_comp(s):
    sc = ''
    for i in reversed(s):
        if i == 'A':
            sc += 'T'
        elif i == 'C':
            sc += 'G'
        elif i == 'G':
            sc += 'C'
        elif i == 'T':
            sc += 'A'
    return sc


if __name__ == "__main__":
    fa = pysam.FastaFile(sys.argv[1])
    res = open(sys.argv[3], "w")
    model = sys.argv[4]
    paths = open(sys.argv[2], 'r')
    # fa = pysam.FastaFile("/Users/troye/Documents/Recycler/bin/sim800/all.fa")
    # res = open("/Users/troye/Documents/Recycler/bin/sim800/neo800.fa", "w")
    # paths = open("/Users/troye/Documents/tmp/seqGraph/r.debug.c.522446.test.txt", 'rb')
    for idx, line in enumerate(paths.readlines()):
        if line.startswith("iter") or line.startswith("self"):
            continue
        fasta = ""
        contigs = line.strip().split('\t')
        for contig in contigs:
            # print(contig)
            if contig[-1] == "+":
                try:
                    fasta += fa.fetch(contig[:-1])
                except Exception as e:
                    name = "_".join(contig[:-1].split('_')[:-1])
                    print("not found: ", name)
                    fasta += fa.fetch(name)
            else:
                try:
                    fasta += rev_comp(fa.fetch(contig[:-1]))
                except Exception as e:
                    name = "_".join(contig[:-1].split('_')[:-1])
                    print("not found: ", name)
                    fasta += rev_comp(fa.fetch(name))
        if model == "0":
            res.write(">res_" + str(idx + 1) + "_" + str(len(fasta))+ "\n")
        else:
            res.write(">"+"".join(contigs)+"\n")
        res.write(fasta + "\n")

    #
    # i = 0
    # for line in paths.readlines():
    #     fasta = ""
    #     if i % 2 != 0:
    #         contigs = line.strip().split('\t')
    #         for contig in contigs:
    #             # print(contig)
    #             if contig[-1]=="+":
    #                 print(i, contig[:-1], contig[-1])
    #                 try:
    #                     fasta += fa.fetch(contig[:-1])
    #                 except Exception as e:
    #                     name = "_".join(contig[:-1].split('_')[:-1])
    #                     print(i, name)
    #                     fasta += fa.fetch(name)
    #             else:
    #                 print(i, contig[:-1], contig[-1])
    #                 try:
    #                     fasta += rev_comp(fa.fetch(contig[:-1]))
    #                 except Exception as e:
    #                     name = "_".join(contig[:-1].split('_')[:-1])
    #                     print(i, name)
    #                     fasta += rev_comp(fa.fetch(name))
    #         res.write(">res"+str(i+1)+"_"+str(len(fasta))+"\n")
    #         res.write(fasta+"\n")
    #     i += 1

    # chx=fa.fetch("EDGE_1298087_length_128_cov_268.000000")
    # print(chx)
    # with pysam.FastxFile("/Users/troye/Documents/Recycler/bin/522446/all.fa") as fin:
    #     chrx= fin.fetch("EDGE_61_length_838_cov_170.571027")
    #     print(chrx)
    # for r in fin:
    #     print(r)
    res.close()
    paths.close()
