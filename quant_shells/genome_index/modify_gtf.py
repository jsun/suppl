import os
import sys
import copy

#
# python modify_gtf.py IWGSC_v1.1_HC_20170706.gff3 > IWGSC_v1.1_HC_20170706.jsun.gff3
# python modify_gtf.py IWGSC_v1.1_LC_20170706.gff3 > IWGSC_v1.1_LC_20170706.jsun.gff3
#
# python modify_gtf.py IWGSC_v1.1_HC_20170706.gtf > IWGSC_v1.1_HC_20170706.jsun.gtf
# python modify_gtf.py IWGSC_v1.1_LC_20170706.gtf > IWGSC_v1.1_LC_20170706.jsun.gtf
#



## The ampping indexes were created from `161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta`.
## But, the coordinates in GTF/GFF (IWGSC v1.1) are designed to `161010_Chinese_Spring_v1.0_pseudomolecules.fasta`.
## This script is used for modifying the coordinates from whle-sequence fasta to parted-sequence fasta by
## using `161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed`.
## The useage of this script is written in the script file.




def load_exdata(filepath):
    exdata = {}
    with open(filepath, 'r') as infh:
        for buf in infh:
            r = buf.replace('\n', '').split('\t')
            if r[4] == '0':
                exdata[r[3]] = int(r[5])
    return exdata




def modify_gff3(filepath, exdata):
    
    with open(filepath, 'r') as infh:
        for buf in infh:
            if buf[0:1] == '#':
                print(buf.replace('\n', ''))
            else:
                r = buf.replace('\n', '').split('\t')
                if int(r[3]) < exdata[r[0]] and int(r[4]) < exdata[r[0]]:
                    r[0] = r[0] + '_part1'
                elif int(r[3]) > exdata[r[0]] and int(r[4]) > exdata[r[0]]:
                    r[3] = str(int(r[3]) - exdata[r[0]])
                    r[4] = str(int(r[4]) - exdata[r[0]])
                    r[0] = r[0] + '_part2'
                elif int(r[3]) < exdata[r[0]] and int(r[4]) > exdata[r[0]]:
                    _chrname = r[0]
                    s = copy.deepcopy(r)
                    r[0] = _chrname + '_part1'
                    r[4] = str(exdata[_chrname])
                    s[0] = _chrname + '_part2'
                    s[3] = '1'
                    s[4] = str(int(s[4]) - exdata[_chrname])
                else:
                    raise ValueError(r)
            
                print('\t'.join(r))






if __name__ == '__main__':
    
    bed_filepath = '161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed'
    exdata = load_exdata(bed_filepath)
    modify_gff3(sys.argv[1], exdata)
    
    
    
