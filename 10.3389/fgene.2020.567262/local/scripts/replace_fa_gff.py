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
    posshift = {}
    for key in genome:
        genomestr[key] = list(genome[key])
        posshift[key] = {}
    return genomestr, posshift

def replace_snps(filename, fasta, posshift):
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
        ps = posshift[cols[0]]
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
        if pos+mutlen in ps:
            sys.stderr.write("DUPLICATE REPLACE! at %d" % (pos+mutlen))
        ps[pos+mutlen] = len(alt) - len(cols[3])
        lendiff = len(alt) - len(seq[pos])
        seq[pos] = alt
        for p in range(pos+1,pos+mutlen):
            sys.stderr.write("REMOVE\t%d, %s\n" % (p, seq[p]))
            lendiff = lendiff - len(seq[p])
            seq[p] = ""
        sys.stderr.write("POS:\t%d, %d, DIFF:%d\n" % (pos, pos+mutlen, lendiff))
        ps[pos+mutlen] = lendiff
        line = fp.readline()
        
    return fasta, posshift

def posshift2posdiff(posshift, fasta):
    chrs = fasta.keys()
    posdiff = {}
    for ch in chrs:
        posdiff[ch] = {}
        ps = posshift[ch]
        chrlen = len(fasta[ch])
        sumdiff = 0
        for i in range(0,chrlen+1):
            if i in ps:
                sumdiff += ps[i]
                sys.stderr.write("SUMDIFF: pos:%s, %d diff:%d\n" % (ch, i, sumdiff))
            posdiff[ch][i] = i + sumdiff
    return posdiff

def replace_gff(gfffile, posdiff, outfile):
    fp = open(gfffile, "r")
    out = open(outfile, "w")
    line = fp.readline()
    while line:
        cols = line.rstrip().split("\t")
        ch = posdiff[cols[0]]
        cols[3] = str(ch[int(cols[3])])
        cols[4] = str(ch[int(cols[4])])
        out.write("\t".join(cols) + "\n")
        line = fp.readline()
    fp.close()
    out.close()


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

def main():
    # 1: in_fasta, 2: vcf, 3: in_gff, 4: out_fasta, 5: out_gff
    fasta, posshift = read_fasta(sys.argv[1])
    repfasta, posshift = replace_snps(sys.argv[2], fasta, posshift)
    fastafp = open(sys.argv[4], "w")
    for f in repfasta:
        fastafp.write(">" + f + "\n")
        fastafp.write(str(formatseq("".join(repfasta[f]))))
    fastafp.close()
    posdiff = posshift2posdiff(posshift, fasta)
    repgff = replace_gff(sys.argv[3], posdiff, sys.argv[5])

if __name__ == '__main__':
    main()
