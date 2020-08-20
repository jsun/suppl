#!/usr/bin/env python
import sys

MINQV = 20
MINDEPTH = int(sys.argv[2])
MAXDEPTH = int(sys.argv[3])
MINALTRATIO = 0.8

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

def main():
    fp = open(sys.argv[1])
    line = fp.readline()
    while line:
        line = line.rstrip()
        if line[0] == "#":
            print line
            line = fp.readline()
            continue
        cols = line.split("\t")
        qual = float(cols[5])
        keyval = kv_decode(cols[7])
        dp4 = keyval["DP4"]
        count = [int(v) for v in dp4.split(",")]
        dp = sum(count)
        if qual < MINQV:
            sys.stderr.write("%s\t%s\tRemove by QV %3.2f\n" % (cols[0], cols[1], qual))
            line = fp.readline()
            continue
        if dp < MINDEPTH or dp > MAXDEPTH:
            sys.stderr.write("%s\t%s\tRemove by depth %s\n" % (cols[0], cols[1], dp))
            line = fp.readline()
            continue
        altratio = float(sum(count[2:4])) / sum(count)
        if altratio < MINALTRATIO:
            sys.stderr.write("%s\t%s\tRemove by ratio %3.2f (%s)\n" % (cols[0], cols[1], altratio, keyval["DP4"]))
            line = fp.readline()
            continue
        if "INDEL" in keyval:
            print line
            line = fp.readline()
            continue
        print line
        line = fp.readline()
        continue

if __name__ == '__main__':
    main()
