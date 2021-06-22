import sys
import os



def seek_gene(ref_fpath, gene_id):
    
    fa2seq = {}
    
    s = False
    
    with open(ref_fpath, 'r') as infh:
        for buf in infh:
            buf = buf.replace('\n', '')
            
            if buf[0] == '>':
                if gene_id in buf:
                    s = True
                else:
                    s = False
                
                
            if s:
                print(buf)
    
    
def seek_genes(ref_fpath):

    genes = ["TraesCS1A02G150900", "TraesCS1A02G400500", "TraesCS2A02G261200",
             "TraesCS2A02G368000", "TraesCS2A02G384300", "TraesCS2A02G543900",
             "TraesCS2A02G588300", "TraesCS3A02G183900", "TraesCS3A02G307900",
             "TraesCS3A02G380400", "TraesCS4A02G177600", "TraesCS5A02G049600",
             "TraesCS5A02G067500", "TraesCS5A02G110900", "TraesCS5A02G137800",
             "TraesCS5A02G349100", "TraesCS5A02G421800", "TraesCS5A02G440800",
             "TraesCS5A02G475600", "TraesCS5A02G478100", "TraesCS6A02G033300",
             "TraesCS6A02G216500", "TraesCS7A02G169100"]
    
    for gene in genes:
        seek_gene(ref_fpath, gene)

    



if __name__ == '__main__':
    
    if len(sys.argv) == 3:
        ref_fpath = sys.argv[1]
        gene_id = sys.argv[2]
    
        seek_gene(ref_fpath, gene_id)
        
    else:
        ref_fpath = sys.argv[1]
        seek_genes(ref_fpath)

