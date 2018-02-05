#!/usr/bin/env python 
import sys

def read_fasta(filename):
    fp = open(filename, "r")
    line = fp.readline()
    genome = {}
    key = ""
    seq = ""
    while line:
        line = line.rstrip()
        if line[0] == ">":
            if seq != "" and key != "":
                genome[key] = seq
                seq = ""
            key = line[1:].split(" ")[0]
        else:
            seq += line
        line = fp.readline()
    if seq != "" and key != "":
        genome[key] = seq
    
    genomestr = {}
    for key in genome:
        genomestr[key] = list(genome[key].lower())
    return genomestr

def add_depth(fastastr, depth, mindepth):
    dfp = open(depth, "r")
    line = dfp.readline()
    count = 0
    while line:
        cols = line.rstrip().split("\t")
        if cols[2] >= mindepth:
            f = fastastr[cols[0]]
	    pos = int(cols[1]) - 1
            f[pos] = f[pos].upper()
        line = dfp.readline()
	count += 1
    return fastastr

def formatseq(seq):
    length = len(seq)
    pos = 0
    skip = 60
    nseq = ""
    while pos < length:
        nseq += seq[pos:pos+skip]
        nseq += "\n"
        pos = pos+skip
    return nseq

def main(fasta, depth, mindepth):
    fastastr = read_fasta(fasta)
    sys.stderr.write("Reading finished\n")
    fastastr = add_depth(fastastr, depth, mindepth)
    for f in fastastr:
        sys.stdout.write(">" + f + "\n")
	#sys.stdout.write(str("".join(fastastr[f])))
        sys.stdout.write(str(formatseq("".join(fastastr[f]))))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
