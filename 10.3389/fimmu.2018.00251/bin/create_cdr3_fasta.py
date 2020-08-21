#!/usr/bin/env python
import os
import sys
import argparse

def create_cdr3aa_fasta(i, fugu_id):
    oprot = i + '.cdr3prot.fa'
    onucl = i + '.cdr3nucl.fa'
    outfhprot = open(oprot, 'w')
    outfhnucl = open(onucl, 'w')
    
    with open(i, 'r') as infh:
        for buf in infh:
            if buf[0:6] == '#BEGIN':
                seq_id    = ''
                seq_op    = '.'
                cdr3_prot = ''
                cdr3_nucl = ''
                est_v     = ''
                est_d     = ''
                est_j     = ''
            if buf[0:2] == 'QN':
                seq_id = fugu_id + '_' + buf[2:].strip()
            if buf[0:2] == 'VN':
                est_v = buf[2:].strip()
            if buf[0:2] == 'DN':
                est_d = buf[2:].strip()
            if buf[0:2] == 'JN':
                est_j = buf[2:].strip()
            if buf[0:2] == 'OP':
                seq_op = buf[3:4]
            if buf[0:2] == 'CA':
                cdr3_nucl = buf[2:].strip()
            if buf[0:7] == '#CDR3AA':
                cdr3_prot = buf[7:].strip()
            if buf[0:4] == '#END':
                if '*' not in cdr3_prot and cdr3_prot != '.':
                    outfhprot.write('>' + seq_id + '\n' + cdr3_prot + '\n')
                    outfhnucl.write('>' + seq_id + '\n' + cdr3_nucl + '\n')
    outfhprot.close()
    outfhnucl.close()
 



def create_vjcdr3nucl(i):
    o1 = i + '.vj_cdr3nucl.txt'
    o2 = i + '.cdr3prot_cdr3nucl.txt'
    
    pydict = {}
    aaseq = {}
    
    N = 0
    with open(i, 'r') as infh:
        for buf in infh:
            if buf[0:6] == '#BEGIN':
                seq_id    = ''
                seq_op    = '.'
                cdr3_prot = ''
                cdr3_nucl = ''
                est_v     = ''
                est_d     = ''
                est_j     = ''
            if buf[0:2] == 'QN':
                seq_id = buf[2:].strip()
            if buf[0:2] == 'VN':
                est_v = buf[2:].strip()
            if buf[0:2] == 'DN':
                est_d = buf[2:].strip()
            if buf[0:2] == 'JN':
                est_j = buf[2:].strip()
            if buf[0:2] == 'OP':
                seq_op = buf[3:4]
            if buf[0:2] == 'CA':
                cdr3_nucl = buf[2:].strip()
            if buf[0:7] == '#CDR3AA':
                cdr3_prot = buf[7:].strip()
            if buf[0:4] == '#END':
                if '*' not in cdr3_prot and len(cdr3_prot) != 0:
                    N += 1
                    vj = est_v + est_j
                    if cdr3_nucl not in pydict:
                        pydict[cdr3_nucl] = {'total': 0, 'combn': {}}
                    pydict[cdr3_nucl]['total'] += 1
                    if vj not in pydict[cdr3_nucl]['combn']:
                        pydict[cdr3_nucl]['combn'][vj] = 0
                    pydict[cdr3_nucl]['combn'][vj] += 1
                    
                    if cdr3_prot not in aaseq:
                        aaseq[cdr3_prot] = {'total': 0, 'ntseq': {}}
                    aaseq[cdr3_prot]['total'] += 1
                    if cdr3_nucl not in aaseq[cdr3_prot]['ntseq']:
                        aaseq[cdr3_prot]['ntseq'][cdr3_nucl] =  0
                    aaseq[cdr3_prot]['ntseq'][cdr3_nucl] += 1
                    
    
    
    with open(o1, 'w') as fh:
        for nucl_seq in pydict.keys():
            txt = nucl_seq + '\t' + str(pydict[nucl_seq]['total']) + '\t'
            txt = txt + str(float(pydict[nucl_seq]['total']) / N) + '\t'
            txt = txt + str(len(pydict[nucl_seq]['combn'])) + '\t'
            for combn in pydict[nucl_seq]['combn'].keys():
                txt = txt + combn + ':' + str(pydict[nucl_seq]['combn'][combn]) + ';'
            fh.write(txt + '\n')
    
    with open(o2, 'w') as fh:
        for prot_seq in aaseq.keys():
            txt = prot_seq + '\t' + str(aaseq[prot_seq]['total']) + '\t'
            txt = txt + str(float(aaseq[prot_seq]['total']) / N) + '\t'
            txt = txt + str(len(aaseq[prot_seq]['ntseq'])) + '\t'
            for nucl_seq in aaseq[prot_seq]['ntseq'].keys():
                txt = txt + nucl_seq + ':' + str(aaseq[prot_seq]['ntseq'][nucl_seq]) + ';'
            fh.write(txt + '\n')

 
    




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Creat CDR3 FASTA File.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-f', '--fugu', required = True)
    args = parser.parse_args()
    create_vjcdr3nucl(args.input)
    create_cdr3aa_fasta(args.input, args.fugu)
    
