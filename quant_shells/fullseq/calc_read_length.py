import os
import sys
import gzip


def calc_len(in_fpath):
    
    out_fpath = in_fpath + '.len.tsv'
    
    i = 0
    with gzip.open(in_fpath, 'rt') as infh, open(out_fpath, 'w') as outfh:
        for fqbuf in infh:
            i = i + 1
            if i == 4:
                outfh.write(str(len(fqbuf) - 1) + '\n')
                i = 0
    




if __name__ == '__main__':
    
    in_fpath = sys.argv[1]
    calc_len(in_fpath)
    
    



