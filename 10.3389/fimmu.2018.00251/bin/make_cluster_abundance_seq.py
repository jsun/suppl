#!/usr/bin/env python
import os
import sys
import argparse
import re
from Bio import SeqIO




def create_tsv(i, o):
    cls2type = {}
    cls2nt = {}
    
    with open(i, 'r') as fhin:
        for buf in fhin:
            r = buf.replace('\r', '').replace('\n', '').split('\t')
            cls2type[r[0]] = r[1]
            if r[0] not in cls2nt:
                cls2nt[r[0]] = {}
            cls2nt[r[0]][r[7]] = 1
    
    with open(o, 'w') as fhout:
        for k in cls2type.keys():
            t = k + '\t' + str(cls2type[k]) + '\t' + str(len(cls2nt[k]))
            fhout.write(t + '\n')
    
    
    cls2fugu2nt = {}
    with open(i, 'r') as fhin:
        for buf in fhin:
            r = buf.replace('\r', '').replace('\n', '').split('\t')
            _cls = r[0]
            _type = r[1]
            _fugu = r[3]
            _seq = r[7]
            
            if _type == '3':
                if _cls not in cls2fugu2nt:
                    cls2fugu2nt[_cls] = {}
                if _fugu not in cls2fugu2nt[_cls]:
                    cls2fugu2nt[_cls][_fugu] = {}
                if _seq not in cls2fugu2nt[_cls][_fugu]:
                    cls2fugu2nt[_cls][_fugu][_seq] = 1
    
    with open(o + '.public', 'w') as fhout:
        for k in cls2fugu2nt:
            _1 = len(cls2fugu2nt[k]['fugu1']) if 'fugu1' in cls2fugu2nt[k] else 0
            _2 = len(cls2fugu2nt[k]['fugu2']) if 'fugu2' in cls2fugu2nt[k] else 0
            _3 = len(cls2fugu2nt[k]['fugu3']) if 'fugu3' in cls2fugu2nt[k] else 0
            _t = int(_1) + int(_2) + int(_3)
            fhout.write(k + '\t' + str(_1) + '\t' + str(_2) + '\t' + str(_3) + '\t' + str(_t) + '\n')
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    create_tsv(args.input, args.output)
    
    
