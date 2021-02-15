import os
import sys
import gzip
import re
import joblib
import glob


RE_POLY_A = r'A{5,}|T{5,}'
MIN_SEQ_LEN = 40
MIN_POLYA_LEN = 10


def seek_polyA_start_pos(seq):
    pos = []
    for m in re.finditer(RE_POLY_A, seq):
        pos.append(m.span())
    return pos


def calc_polyA_length(pos):
    polyA_len = []
    if len(pos) > 0:
        for _pos in pos:
            polyA_len.append(_pos[1] - _pos[0])
    return polyA_len


def trim_polyA(input_file_path, output_file_path):
    
    with gzip.open(input_file_path, 'rt') as infh, gzip.open(output_file_path, 'wt') as outfh:
        i = 0
        pos = -1
        pos_list = []
        seq = ''
        seq_len = 0
        for buf in infh:
            i = i + 1
            if i == 1:
                seq = seq + buf
            elif i == 3:
                seq = seq + buf
            elif i == 2:
                pos_list = seek_polyA_start_pos(buf)
                polyA_len = calc_polyA_length(pos_list)
                pos = len(buf[: -1])
                if len(pos_list) > 0:
                    if max(polyA_len) > MIN_POLYA_LEN or sum(polyA_len) > (len(buf[:-1]) / 3):
                        pos = pos_list[0][0]
                
                seq = seq + buf[:pos] + '\n'
                seq_len = len(buf[:pos])
            elif i == 4:
                seq = seq + buf[:pos] + '\n'
                if seq_len >= MIN_SEQ_LEN:
                    outfh.write(seq)
                
                i = 0
                n_trim = 0
                seq = ''
                seq_len = 0
                pos = -1
                pos_list = []
 



def process_fastq(input_dpath=None, output_dpath=None, input_pattern=None, output_pattern=None, n_cpu=-1):
    
    if os.path.isfile(input_dpath):
        trim_polyA(input_dpath, output_dpath)
    
    else:
        
        fq_input_fpath = sorted(glob.glob(os.path.join(input_dpath, '*' + input_pattern)))
        fq_output_fpath = []
        for f in fq_input_fpath:
            fq_output_fpath.append(f.replace(input_pattern, output_pattern))
        
        r = joblib.Parallel(n_jobs = n_cpu)([joblib.delayed(trim_polyA)(in_fpath, out_fpath) \
                                                for in_fpath, out_fpath in zip(fq_input_fpath, fq_output_fpath)])
        
        
    


       
if __name__ == '__main__':
    
    input_dpath = sys.argv[1]
    output_dpath = sys.argv[2]
    if len(sys.argv) > 3:
        input_pattern = sys.argv[3]
        output_pattern = sys.argv[4]
        n_cpu = int(sys.argv[5])
    else:
        input_pattern = None
        output_pattern = None
        n_cpu = 1
        
    process_fastq(input_dpath, output_dpath, input_pattern, output_pattern, n_cpu)





