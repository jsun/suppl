import sys
import os
import gzip
import random

SEQLEN = list(range(40, 107))
SEQPROB = [0.011731536, 0.008515791, 0.008423760, 0.008973900, 0.008520556, 0.010943239,
           0.008820157, 0.009040789, 0.009108830, 0.009083783, 0.010427375, 0.009211721,
           0.009667194, 0.009955498, 0.010696711, 0.011133154, 0.010929137, 0.009652486,
           0.008298355, 0.007942544, 0.007285223, 0.006761953, 0.006790285, 0.006547897,
           0.006549830, 0.005963953, 0.006098910, 0.006323993, 0.006286632, 0.005951110,
           0.005713492, 0.006544578, 0.005962626, 0.006027030, 0.006443726, 0.006348507,
           0.006453347, 0.006118867, 0.006458880, 0.006996158, 0.006798700, 0.007557760,
           0.007489487, 0.007104454, 0.006941834, 0.007121398, 0.007842928, 0.008183169,
           0.007694717, 0.008575665, 0.009022032, 0.009653639, 0.008868935, 0.009730585,
           0.010179781, 0.009787761, 0.008274973, 0.008572619, 0.008435650, 0.008118566,
           0.008618231, 0.008394715, 0.009092257, 0.011835133, 0.026041551, 0.113266789,
           0.338093158]


def shortage_fq(fq_fpath, fq_out):
    random_seed = abs(hash(os.path.basename(fq_fpath))) % (10 ** 8)
    random.seed(random_seed)
    
    
    with gzip.open(fq_fpath, 'rt') as infh, gzip.open(fq_out, 'wt') as outfh:
        i = 0
        fq_dat = []
        for fqbuff in infh:
            i = i + 1
            fqbuff = fqbuff.replace('\n', '')
            fq_dat.append(fqbuff)
            
            if i == 4:
                seqlen = random.choices(SEQLEN, weights=SEQPROB, k=1)[0]
                fq_dat[1] = fq_dat[1][0:seqlen]
                fq_dat[3] = fq_dat[3][0:seqlen]
                
                outfh.write('\n'.join(fq_dat) + '\n')
                i = 0
                fq_dat = []
    




if __name__ == '__main__':
    fq_fpath = sys.argv[1]
    fq_out = fq_fpath.replace('fastq.gz', '') + 'short.fastq.gz'
    
    shortage_fq(fq_fpath, fq_out)






