#!/usr/bin/env python
import sys

def kv_decode(form):
    cols = form.split(";")
    keyval = {}
    for c in cols:
        kv = c.split("=")
        if len(kv) == 1:
            keyval[kv[0]] = "T"
        else:
            keyval[kv[0]] = kv[1]
    return keyval

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
        genomestr[key] = list(genome[key])
    return genomestr

def replace_snps(filename, fasta):
    fp = open(filename, "r")
    line = fp.readline()
    while line:
        line = line.rstrip()
        if line[0] == "#":
            line = fp.readline()
            continue
        cols = line.split("\t")
        keyval = kv_decode(cols[7])
        seq = fasta[cols[0]]
        pos = int(cols[1])-1
        mutlen = len(cols[3])
        alt = cols[4]
        qual = float(cols[5])
        sys.stderr.write("RAW\t%s\n" % line)
        if alt.find(",") >= 0:
            alts = alt.split(",")
            if cols[3] in alts:
                sys.stderr.write("Warning\tSkip due to multiple replacement candidates\n")
                line = fp.readline()
                continue
            alt = alts[0]
        if "INDEL" in keyval and qual < 50:
            sys.stderr.write("Warning\tSkip due to low quality INDEL\n")
            line = fp.readline()
            continue
        if seq[pos] == "":
            sys.stderr.write("Warning\tSkip this mutation because of close to INDEL site\n")
            line = fp.readline()
            continue
        if "".join(seq[pos:(pos+mutlen)]) != cols[3]:
            sys.stderr.write("Warning\t%s is not %s (%s) at %s, %s\n" % 
                             (cols[3], "".join(seq[pos:(pos+mutlen)]), 
                              "".join(seq[(pos-1):(pos+mutlen+1)]), 
                              cols[0], cols[1]))
        sys.stderr.write("REPLACE\t%d\t%s\t%s\t%s\n" % (pos, seq[pos], cols[3], alt))
        seq[pos] = alt
        for p in range(pos+1,pos+mutlen):
            sys.stderr.write("REMOVE\t%d\n" % p)
            seq[p] = ""
        line = fp.readline()
        
    return fasta

def main():
    fasta = read_fasta(sys.argv[1])
    repfasta = replace_snps(sys.argv[2], fasta)
    for f in repfasta:
        print ">" + f
        print "".join(repfasta[f])

if __name__ == '__main__':
    main()
