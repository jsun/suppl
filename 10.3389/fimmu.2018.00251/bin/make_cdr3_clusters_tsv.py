#!/usr/bin/env python
import os
import sys
import argparse
import re
from Bio import SeqIO


def create_tsv(f, c, o):
    # read FASTA file
    fastafh = open(f, 'rU')
    seqid2seq = {}
    for record in SeqIO.parse(fastafh, "fasta"):
        seqid2seq[record.description] = str(record.seq)
    fastafh.close()
    
    
    # read cluster file
    mptn = re.compile(r'>(.+)\.\.\.')
    clusterid2seqidlist     = {}
    clusterid2repseqid      = {}
    clusterid2seqsize_fugu1 = {}
    clusterid2seqsize_fugu2 = {}
    clusterid2seqsize_fugu3 = {}
    cluster_id = None
    with open(c, 'r') as clstrfh:
        for buf in clstrfh:
            if buf[0:1] == '>':
                cluster_id = buf[1:]
                cluster_id = cluster_id.replace('\n', '')
                cluster_id = cluster_id.replace(' ', '_')
                # init
                clusterid2seqidlist[cluster_id]     = []
                clusterid2repseqid[cluster_id]      = None
                clusterid2seqsize_fugu1[cluster_id] = 0
                clusterid2seqsize_fugu2[cluster_id] = 0
                clusterid2seqsize_fugu3[cluster_id] = 0
            else:
                m = mptn.search(buf)
                seqid = m.group(1)
                if 'fugu1' in buf:
                    clusterid2seqsize_fugu1[cluster_id] += 1
                if 'fugu2' in buf:
                    clusterid2seqsize_fugu2[cluster_id] += 1
                if 'fugu3' in buf:
                    clusterid2seqsize_fugu3[cluster_id] += 1
                if '*' in buf:
                    clusterid2seqidlist[cluster_id].append('*' + seqid)
                    clusterid2repseqid[cluster_id] = seqid
                else:
                    clusterid2seqidlist[cluster_id].append(seqid)
    
    # print out tsv
    with open(o, 'w') as outfh:
        outfh.write('ClusterID\tRepresentSeq\tRepresentSeqLen\tFugu1Count\tFugu2Count\tFugu3Count\tTotalCount\tSeqID')
        for cls_id in sorted(clusterid2repseqid.iterkeys()):
            arr = [cls_id,
                   seqid2seq[clusterid2repseqid[cls_id]],
                   str(len(seqid2seq[clusterid2repseqid[cls_id]])),
                   str(clusterid2seqsize_fugu1[cls_id]),
                   str(clusterid2seqsize_fugu2[cls_id]),
                   str(clusterid2seqsize_fugu3[cls_id]),
                   str(clusterid2seqsize_fugu1[cls_id] + clusterid2seqsize_fugu2[cls_id] + clusterid2seqsize_fugu3[cls_id]),
                   ';'.join(clusterid2seqidlist[cls_id])]
            outfh.write('\t'.join(arr) + '\n')
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Create CDR3 capture-recapture dataset.')
    parser.add_argument('-f', '--fasta', required = True)
    parser.add_argument('-c', '--clstr', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    create_tsv(args.fasta, args.clstr, args.output)
    
    
