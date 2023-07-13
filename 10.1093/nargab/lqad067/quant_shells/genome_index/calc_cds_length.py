from pandabox import seqUtils



def calc_length(fpath):
    
    fx = seqUtils.FASTX()
    
    for fx_ in fx.parse_fasta(fpath):
        seq_id = fx_['id'].split(' ')[0]
        print('{}\t{}'.format(seq_id, len(fx_['seq'])))


if __name__ == '__main__':
    apath = 'chrA.cds.fa'
    bpath = 'chrB.cds.fa'
    dpath = 'chrD.cds.fa'

    calc_length(apath)
    calc_length(bpath)
    calc_length(dpath)
    
    

