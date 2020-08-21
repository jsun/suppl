#!/usr/bin/env python
import os
import sys
import argparse


phash = {'nFVH1': 'v1', 'nFVH2': 'v2', 'nFVH3': 'v3',
         'nVhCm1': 'cm', 'nVhCt1': 'ct'}


def read_cutadaptlog(logf, pdict):
    with open(logf, 'r') as infh:
        for buf in infh:
            rd = buf.split('\t')
            if len(rd) > 7:
                rd2 = rd[0].split(' ')
                seqid = rd2[0]
                primer = rd[7]
                if seqid not in pdict:
                    pdict[seqid] = {'V': None, 'C': None}
                if 'nFVH' in primer:
                    pdict[seqid]['V'] = phash[primer]
                elif 'nVhC' in primer:
                    pdict[seqid]['C'] = phash[primer]
    return pdict
    

def printout_fq(fq, pdict):
    fqout_v1cm = fq[:-3] + '.v1.cm.fq'
    fqout_v2cm = fq[:-3] + '.v2.cm.fq'
    fqout_v3cm = fq[:-3] + '.v3.cm.fq'
    fqout_v1ct = fq[:-3] + '.v1.ct.fq'
    fqout_v2ct = fq[:-3] + '.v2.ct.fq'
    fqout_v3ct = fq[:-3] + '.v3.ct.fq'
    
    fqoutfh = {}
    fqoutfh['v1cm'] = open(fqout_v1cm, 'w')
    fqoutfh['v2cm'] = open(fqout_v2cm, 'w')
    fqoutfh['v3cm'] = open(fqout_v3cm, 'w')
    fqoutfh['v1ct'] = open(fqout_v1ct, 'w')
    fqoutfh['v2ct'] = open(fqout_v2ct, 'w')
    fqoutfh['v3ct'] = open(fqout_v3ct, 'w')
    
    with open(fq, 'r') as fqin:
        cnt = 0
        seqid  = None
        seqdat = None
        for buf in fqin:
            cnt += 1
            if cnt % 4 == 1:
                if seqdat is not None:
                    if seqid in pdict:
                        if pdict[seqid]['V'] is not None and pdict[seqid]['C'] is not None:
                            fqoutfh[pdict[seqid]['V'] + pdict[seqid]['C']].write(seqdat)
                seqid  = buf.split(' ')[0][1:]
                seqdat = buf
            else:
                seqdat = seqdat + buf
    
    fqoutfh['v1cm'].close()
    fqoutfh['v2cm'].close()
    fqoutfh['v3cm'].close()
    fqoutfh['v1ct'].close()
    fqoutfh['v2ct'].close()
    fqoutfh['v3ct'].close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Split FASTQ file according to primers.')
    parser.add_argument('-a', '--fq1', required = True)
    parser.add_argument('-b', '--fq2', required = True)
    parser.add_argument('-c', '--log1', required = True)
    parser.add_argument('-d', '--log2', required = True)
    args = parser.parse_args()
    
    pdict = {}
    pdict = read_cutadaptlog(args.log1, pdict)
    pdict = read_cutadaptlog(args.log2, pdict)
    printout_fq(args.fq1, pdict)
    printout_fq(args.fq2, pdict)
    
    
    
