#!/usr/bin/env python
import os
import sys
import argparse
import re
from Bio import SeqIO



def create_tsv(p, c, o):
    pfiles = p.split(",")
    
    mptn = re.compile(r'>(fugu[1-3])_(.+)\.\.\.')
    
    class_type = {}
    seq_dict = {}
    
    with open(c, 'r') as clstrfh:
        cluster_id = None
        fugu_id = None
        seq_id = None
        for buf in clstrfh:
            if buf[0:1] == '>':
                cluster_id = buf[1:]
                cluster_id = cluster_id.replace('\n', '')
                cluster_id = cluster_id.replace(' ', '_')
            else:
                m = mptn.search(buf)
                fugu_id = m.group(1)
                seq_id = m.group(2)
                seq_dict[seq_id] = {'F': fugu_id, 'C': cluster_id}
                if cluster_id not in class_type:
                    class_type[cluster_id] = {'N': 0}
                class_type[cluster_id][fugu_id] = 1
                class_type[cluster_id]['N'] += 1
    
    outfh = open(o, 'w')
    
    for pfile in pfiles:
        with open(pfile, 'r') as pfh:
            for buf in pfh:
                buf = buf.replace('\n', '')
                if buf[0:6] == '#BEGIN':
                    seq_id = None
                    vdel = None
                    jdel = None
                    vjins = None
                    cdr3aa = None
                if buf[0:2] == 'QN':
                    seq_id = buf[3:]
                if buf[0:2] == 'VD':
                    vdel = buf[3:]
                    if vdel == '.':
                        vdel = ''
                if buf[0:2] == 'JD':
                    jdel = buf[3:]
                    if jdel == '.':
                        jdel = ''
                if buf[0:2] == 'VJ':
                    vjins = buf[3:]
                    if vjins == '.':
                        vjins = ''
                if buf[0:7] == '#CDR3AA':
                    cdr3aa = buf[8:]
                if buf[0:4] == '#END':
                    if seq_id in seq_dict:
                        txt = seq_dict[seq_id]['C'] + '\t' + str(len(class_type[seq_dict[seq_id]['C']]) - 1) + '\t'
                        txt = txt + str(class_type[seq_dict[seq_id]['C']]['N']) + '\t'
                        txt = txt + seq_dict[seq_id]['F'] + '\t'# + seq_id + '\t' 
                        txt = txt + cdr3aa + '\t' + vdel + '\t' + jdel + '\t' + vjins + '\n'
                        outfh.write(txt)
    outfh.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Create CDR3 capture-recapture dataset.')
    parser.add_argument('-p', '--pydair', required = True)
    parser.add_argument('-c', '--clstr', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    create_tsv(args.pydair, args.clstr, args.output)
    
    
