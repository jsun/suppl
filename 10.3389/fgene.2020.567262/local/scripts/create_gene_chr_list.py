import os
import sys
import re
import argparse

#
# ptyhon $this --gff camara/genome.modified.gff 
#
    
def proc_main(gff_path):
    gene2chr = {}
    with open(gff_path, 'r') as fh:
        for buf in fh:
            r = buf.split('\t')
            if r[2] == 'gene':
                gene2chr[r[8][3:14]] = r[0]
    
    for k, v in gene2chr.items():
        print(k + '\t' + v)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'GO data generation.')
    parser.add_argument('-g', '--gff', required = True)
    args = parser.parse_args()
    proc_main(args.gff)
    
