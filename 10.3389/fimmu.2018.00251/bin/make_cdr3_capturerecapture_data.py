#!/usr/bin/env python
import os
import sys
import argparse
import re


def create_dataset(i, o):
    mptn = re.compile(r'>(.+)\.\.\.')
    cluster_id = None
    outfh = open(o, 'w')
    with open(i, 'r') as infh:
        for buf in infh:
            if buf[0:1] == '>':
                cluster_id = buf[1:]
                cluster_id = cluster_id.replace('\n', '')
                cluster_id = cluster_id.replace(' ', '_')
            else:
                m = mptn.search(buf)
                s = m.group(1).split('_', 1)
                outfh.write(cluster_id + '\t' + s[0] + '\t' + m.group(1) + '\n')
    outfh.close()







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Create CDR3 capture-recapture dataset.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    create_dataset_all(args.input, args.output)
    
    
