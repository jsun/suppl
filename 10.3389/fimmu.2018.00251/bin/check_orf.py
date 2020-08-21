#!/usr/bin/env python
import os
import sys
import argparse

def parse_pydair_simple(i):
    cnt_productive = 0
    cnt_nonproductive = 0
    cnt_cdr3len0 = 0
    cnt_cdr3orf = 0
    cnt_cdr3stop = 0
    with open(i, 'r') as fh:
        for buf in fh:
            if buf[0:2] == 'OP':
                orf_codes = buf[3:4]
                if orf_codes in ['1', '2', '3']:
                    has_orf = True
                    cnt_productive += 1
                else:
                    cnt_nonproductive += 1
            if buf[0:5] == '#CDR3':
                cdr3aa = buf[8:].rstrip()
                if cdr3aa == '.':
                    cnt_cdr3len0 += 1
                elif '*' in cdr3aa:
                    cnt_cdr3stop += 1
                else:
                    cnt_cdr3orf += 1
                
    print [cnt_productive, cnt_nonproductive,
           cnt_cdr3orf, cnt_cdr3stop, cnt_cdr3len0]
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Check the ORF of CDR3 sequences.')
    parser.add_argument('-i', '--input', required = True)
    args = parser.parse_args()
    print args.input
    parse_pydair_simple(args.input)



