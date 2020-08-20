#!/usr/bin/env python
import sys

def main():
    fp = open(sys.argv[1], "r")
    vcfs = []
    for line in fp:
        if line[0] == "#":
            print line.rstrip()
            continue
        cols = line.split("\t")
        vcfs.append([cols[0], int(cols[1]), line.rstrip()])
    vcfs = sorted(vcfs, key=lambda x:(x[0],x[1]))
    for v in vcfs:
        print v[2]

if __name__ == '__main__':
    main()
