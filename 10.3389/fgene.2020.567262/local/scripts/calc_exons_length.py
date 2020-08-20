import os
import sys
import re
import argparse

#
# ptyhon $this --gffa camara/genome.modified.gff --gffr crivularis/genome.modified.gff --output camara.gene.length
#

def proc_main(gff_file):
    
    idptn = re.compile("ID=(CARHR[0-9]+)\..+")
    
    gene_length = {}
    with open(gff_file, 'r') as gfffh:
        for buf in gfffh:
            buf_records = buf.split('\t')
            if buf_records[2] == 'exon':
                mat = idptn.match(buf_records[8])
                if mat:
                    gene_name = mat.group(1)
                    if gene_name not in gene_length:
                        gene_length[gene_name] = []
                    gene_length[gene_name].append([int(buf_records[3]), int(buf_records[4])])
    
    
    gene_length_max = {}
    gene_length_min = {}
    for gene_name in gene_length.keys():
        for exon_range in gene_length[gene_name]:
            # define
            if gene_name not in gene_length_max:
                gene_length_max[gene_name] = max(exon_range)
            if gene_name not in gene_length_min:
                gene_length_min[gene_name] = min(exon_range)
            # updaet
            if max(exon_range) > gene_length_max[gene_name]:
                gene_length_max[gene_name] = max(exon_range)
            if min(exon_range) < gene_length_min[gene_name]:
                gene_length_min[gene_name] = min(exon_range)
    
    gene_length_bits = {}
    for gene_name in gene_length.keys():
        if gene_name not in gene_length_bits:
            gene_length_bits[gene_name] = [0] * (gene_length_max[gene_name] - gene_length_min[gene_name] + 1)
        for exon_range in gene_length[gene_name]:
            for i in range(min(exon_range), max(exon_range) + 1):
                pos_i = i - gene_length_min[gene_name]
                gene_length_bits[gene_name][pos_i] = 1
   
    gene_length_final = {}
    for gene_name in gene_length.keys():
        gene_length_final[gene_name] = sum(gene_length_bits[gene_name])
    
    return gene_length_final

    
    
    
def write_data(gene_length_ama, gene_length_riv, output):
    with open(output, 'w') as fh:
        for gene_name in gene_length_ama.keys():
            fh.write(gene_name + '\t' + str(gene_length_ama[gene_name]) + '\t' + str(gene_length_riv[gene_name]) + '\n')
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Calculate the length of overlaps of exons for each gene.')
    parser.add_argument('-a', '--gffa', required = True)
    parser.add_argument('-r', '--gffr', required = True)
    parser.add_argument('-g', '--output', required = True)
    args = parser.parse_args()
    gene_length_ama = proc_main(args.gffa)
    gene_length_riv = proc_main(args.gffr)
    write_data(gene_length_ama, gene_length_riv, args.output)
    
    
    
    
    
    
