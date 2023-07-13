import sys
import os


def make_subgenome_fasta(file_path, out_path_A, out_path_B, out_path_D, out_path_unknown):
    
    isChr = None
    
    with open(file_path, 'r') as infh:
        with open(out_path_A, 'w') as outA, open(out_path_B, 'w') as outB, open(out_path_D, 'w') as outD, open(out_path_unknown, 'w') as outU:
            for buf in infh:
                if buf[0:4] == '>chr':
                    isChr = buf[5:6]
                
                if isChr == 'A':
                    outA.write(buf)
                elif isChr == 'B':
                    outB.write(buf)
                elif isChr == 'D':
                    outD.write(buf)
                else:
                    outU.write(buf)
            



if __name__ == '__main__':
    make_subgenome_fasta(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    
    

