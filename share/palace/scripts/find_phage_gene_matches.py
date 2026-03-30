#!/usr/bin/env python
import argparse, glob, os, subprocess, re, shutil, sys
from subprocess import CalledProcessError

def parse_user_input():
    parser = argparse.ArgumentParser(description='find_phage_gene_matches supporting BLAST+, MMseqs2 or DIAMOND')
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('-f','--fasta', help='input file - contigs in fasta format', type=str)
    g.add_argument('-db', '--dbpath', help='path of existing database for the contigs', type=str)
    parser.add_argument('-o','--output', help='output directory', required=True, type=str)
    parser.add_argument('-g','--genes_dir', help='dir with gene nt fasta files', type=str)
    parser.add_argument('-p','--proteins_dir', help='dir with protein aa fasta files', type=str)
    parser.add_argument('--engine', choices=['blast', 'mmseqs', 'diamond'], default='blast', help='Search engine')
    parser.add_argument('--bin_path', help='bin path', default="", type=str)
    parser.add_argument('-n', '--nthreads', help='threads', default=1, type=int)
    parser.add_argument('-t', '--thresh', help='threshold', default=0.75, type=float)
    parser.add_argument('-c', '--clean', help='remove intermediate files', action='store_true')
    return parser.parse_args()

def run_cmd(command):
    print(f"Running command: {command}")
    try:
        subprocess.check_call(command, shell=True)
    except CalledProcessError as e:
        print(f"Error executing command.")
        raise e

def create_db(infile, engine, bin_path, outdir):
    basename = os.path.basename(infile)
    if engine == 'blast':
        dbpath = os.path.join(outdir, basename + ".blastdb")
        cmd = f"{os.path.join(bin_path, 'makeblastdb')} -in {infile} -dbtype nucl -out {dbpath}"
        run_cmd(cmd)
    elif engine == 'diamond':
        return infile
    else:
        dbpath = os.path.join(outdir, basename + ".mmseqsdb")
        cmd = f"{os.path.join(bin_path, 'mmseqs')} createdb {infile} {dbpath}"
        run_cmd(cmd)
    return dbpath

def mmseqs_search_logic(query_fasta, target_db, numthreads, outdir, bin_path, search_type):
    base = os.path.basename(query_fasta)
    tmp_dir = os.path.join(outdir, f"tmp_{base}")
    if not os.path.exists(tmp_dir): os.makedirs(tmp_dir)

    output_txt = os.path.join(outdir, f"{base}_mmseqs.out")
    q_db = os.path.join(tmp_dir, "query.db")
    run_cmd(f"{os.path.join(bin_path, 'mmseqs')} createdb {query_fasta} {q_db}")

    res_db = os.path.join(tmp_dir, "res.db")
    s_cmd = (f"{os.path.join(bin_path, 'mmseqs')} search {q_db} {target_db} {res_db} {tmp_dir} "
             f"--threads {numthreads} --search-type {search_type} -s 4.0 --min-seq-id 0.7")
    run_cmd(s_cmd)

    fmt = "query,target,alnlen,pident,qlen,tlen,evalue"
    c_cmd = (f"{os.path.join(bin_path, 'mmseqs')} convertalis {q_db} {target_db} {res_db} {output_txt} "
             f"--format-output {fmt}")
    run_cmd(c_cmd)

    return output_txt

def search_protein(pf, dbpath, engine, bin_path, numthreads, outdir):
    basename = os.path.basename(pf)
    if engine == 'blast':
        outputpath = os.path.join(outdir, basename + "_blast.out")
        cmd = (f"{os.path.join(bin_path, 'tblastn')} -db {dbpath} -db_gencode 11 -query {pf} "
               f"-out {outputpath} -num_threads {numthreads} "
               f'-outfmt "6 qseqid sseqid length pident qlen slen evalue"')
        run_cmd(cmd)
        return outputpath
    elif engine == 'diamond':
        outputpath = os.path.join(outdir, basename + "_diamond.out")
        dbpath_dmnd = os.path.join(outdir, basename + ".dmnd")
        
        run_cmd(f"{os.path.join(bin_path, 'diamond')} makedb --in {pf} -d {dbpath_dmnd} --quiet")
        
        cmd = (f"{os.path.join(bin_path, 'diamond')} blastx -d {dbpath_dmnd} -q {dbpath} "
               f"-o {outputpath} -p {numthreads} "
               f'--outfmt 6 sseqid qseqid length pident slen qlen evalue')
        run_cmd(cmd)
        
        if os.path.exists(dbpath_dmnd): os.remove(dbpath_dmnd)
        return outputpath
    else:
        return mmseqs_search_logic(pf, dbpath, numthreads, outdir, bin_path, 2)

def search_gene(gf, dbpath, engine, bin_path, numthreads, outdir):
    basename = os.path.basename(gf)
    if engine == 'diamond':
        print(f"Warning: DIAMOND does not support nt vs nt search. Skipping {basename}...")
        return None
        
    if engine == 'blast':
        outputpath = os.path.join(outdir, basename + "_blast.out")
        cmd = (f"{os.path.join(bin_path, 'blastn')} -task megablast -db {dbpath} -query {gf} "
               f"-out {outputpath} -num_threads {numthreads} "
               f'-outfmt "6 qseqid sseqid length pident qlen slen evalue"')
        run_cmd(cmd)
        return outputpath
    else:
        return mmseqs_search_logic(gf, dbpath, numthreads, outdir, bin_path, 3)

def get_hits(blastfile, hit_contigs_set, thresh=0.75, is_protein=False, engine='blast'):
    if not blastfile or not os.path.exists(blastfile): return
    with open(blastfile) as f:
        for line in f:
            splt = line.strip().split('\t')
            if len(splt) < 7: continue

            contig = re.split(':|;', splt[1])[0]
            percentid = float(splt[3])
            matchlen = int(splt[2])
            genelen = int(splt[4])

            if engine == 'mmseqs' and is_protein:
                matchlen = matchlen / 3.0

            coverage = float(matchlen) / float(genelen)

            if percentid > thresh * 100 and coverage > thresh:
                hit_contigs_set[contig] = hit_contigs_set.get(contig, 0) + 1

def main():
    args = parse_user_input()
    
    if args.engine == 'diamond' and not args.fasta:
        print("Error: DIAMOND engine requires the original contigs fasta (-f) to act as the query.")
        sys.exit(1)
        
    if not os.path.exists(args.output): os.makedirs(args.output)
    
    dbpath = args.dbpath if args.dbpath else create_db(args.fasta, args.engine, args.bin_path, args.output)
    results_info = []

    if args.genes_dir and os.path.exists(args.genes_dir):
        for fname in os.listdir(args.genes_dir):
            r_file = search_gene(os.path.join(args.genes_dir, fname), dbpath, args.engine, args.bin_path, args.nthreads, args.output)
            if r_file: results_info.append((r_file, False))

    if args.proteins_dir and os.path.exists(args.proteins_dir):
        for fname in os.listdir(args.proteins_dir):
            r_file = search_protein(os.path.join(args.proteins_dir, fname), dbpath, args.engine, args.bin_path, args.nthreads, args.output)
            if r_file: results_info.append((r_file, True))

    hit_contigs_set = {}
    for r_file, is_prot in results_info:
        get_hits(r_file, hit_contigs_set, args.thresh, is_protein=is_prot, engine=args.engine)

    with open(os.path.join(args.output, "hit_seqs.out"), 'w') as f:
        for k, v in hit_contigs_set.items(): f.write(f"{k}\t{v}\n")

    if args.clean:
        for d in glob.glob(os.path.join(args.output, "tmp_*")): shutil.rmtree(d)
        if args.fasta and args.engine != 'diamond':
            for f in glob.glob(dbpath + "*"): os.remove(f)

if __name__ == '__main__':
    main()
