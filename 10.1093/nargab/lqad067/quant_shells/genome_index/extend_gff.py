import os
import sys
import re
import math

N_EXTEND = None


def load_gff(input_gff):
    
    gff = []
    with open(input_gff, 'r') as infh:
        for buf in infh:
            if buf[0] == '#':
                continue
            
            buf_records = buf.split('\t')
            
            if buf_records[2] == 'gene':
                _gff = {
                    'gene'  : buf_records[8].split(';')[0].replace('ID=', '').replace('\n', ''),
                    'left'  : int(buf_records[3]),
                    'right' : int(buf_records[4]),
                    'strand': buf_records[6],
                }
                
               
                gff.append(_gff)
    
    return gff
    


def calc_interval(gff, i):
    
    n_interval = 0
    
    if i == 0:
        n_interval = gff[i]['left'] - N_EXTEND

    elif i == len(gff) and gff[i]['strand'] == '+':
        n_interval = N_EXTEND

    else:
        if gff[i]['strand'] == '+':
            for j in range(i + 1, len(gff)):
                if gff[i]['right'] < gff[j]['left']:
                    n_interval = gff[j]['left'] - gff[i]['right']
                    break
        else:
            for j in reversed(range(i)):
                if gff[j]['right'] < gff[i]['left']:
                    n_interval = gff[i]['left'] - gff[j]['right']
                    break
    
    if n_interval < 0:
        n_interval = gff[i]['left'] - 1
    
    if n_interval > N_EXTEND:
        n_interval = N_EXTEND
    
    return n_interval
    
    



def extend_genes(gff):
    
    
    for i in range(len(gff)):
        
        n_extend = calc_interval(gff, i)
        
        if gff[i]['strand'] == '+':
            gff[i]['edited_left'] = gff[i]['left']
            gff[i]['edited_right'] = gff[i]['right'] + n_extend
            
        elif gff[i]['strand'] == '-':
            gff[i]['edited_left'] = gff[i]['left'] - n_extend
            gff[i]['edited_right'] = gff[i]['right']

        else:
            print(buf)
            raise ValueError('something wrong ......')
    
    
    return gff




def __gff_iter(input_gff):
    
    gff_records = None
    
    with open(input_gff, 'r') as infh:
        for buf in infh:
            if buf[0] == '#':
                continue
            
            buf = buf.replace('\n', '')
            buf_records = buf.split('\t')
            if buf_records[2] == 'gene':
                if gff_records is not None:
                    yield gff_records
                gff_records = []
            
            gff_records.append(buf)
    
    yield gff_records




def get_idx(gff_records, ftype):
    idx = []
    for i, gff_record in enumerate(gff_records):
        if gff_record.split('\t')[2] == ftype:
            idx.append(i)
    
    return idx





def build_dict(gff):
    gff_dict = {}
    for gff_record in gff:
        gff_dict[gff_record['gene']] = gff_record

    return gff_dict






def modify_gff(input_gff, gff_dict):
    
    for gff_records in __gff_iter(input_gff):
        
        # gene id
        gene_records = gff_records[0].split('\t')
        gene_id = gene_records[8].split(';')[0].replace('ID=', '')
        
        # update gene record
        gene_records[3] = str(gff_dict[gene_id]['edited_left'])
        gene_records[4] = str(gff_dict[gene_id]['edited_right'])
        gff_records[0] = '\t'.join(gene_records)
        
        
        # find the indexs of exons
        idx_exons = get_idx(gff_records, 'exon')
        
        # update exon record
        if gff_dict[gene_id]['strand'] == '+':
            
            # get the last exons and add a new exon after that
            max_idx_exons = max(idx_exons)
        
            # get the last exons and copied information
            r_last = gff_records[max_idx_exons]
            r_extend = r_last.split('\t')
        
            r_extend[3] = str(int(r_last.split('\t')[4]) + 1)
            r_extend[4] = str(gff_dict[gene_id]['edited_right'])
            r_extend[8] = re.sub('_exon_[0-9]+;', '_exon_ext;', r_extend[8])
            gff_records.insert(max_idx_exons + 1, '\t'.join(r_extend))
        
        elif gff_dict[gene_id]['strand'] == '-':
            
            # get the first exons and add a new exon before that
            min_idx_exons = min(idx_exons)
        
            # get the last exons and copied information
            r_first = gff_records[min_idx_exons]
            r_extend = r_first.split('\t')
            
            r_extend[3] = str(gff_dict[gene_id]['edited_left'])
            r_extend[4] = str(int(r_first.split('\t')[3]) - 1)
            r_extend[8] = re.sub('_exon_[0-9]+;', '_exon_ext;', r_extend[8])
            gff_records.insert(min_idx_exons, '\t'.join(r_extend))
            
        else:
            raise ValueError('Unknown strand ....')
        
        
        print('\n'.join(gff_records))
    


if __name__ == '__main__':
    
    input_gff = sys.argv[1]
    N_EXTEND = int(sys.argv[2])
    
    gff = load_gff(input_gff)
    gff = extend_genes(gff)
    gff_dict = build_dict(gff)
    modify_gff(input_gff, gff_dict)


