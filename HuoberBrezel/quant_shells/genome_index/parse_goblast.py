import os
import sys
import copy

#
# python $this blastresult.tsv golist.txt
#



# TraesCS1A02G064900.2    1428    LOC_Os02g36440.1        1665    80.192  313     60      2       5       316     152     463     5.09e-60        233
# TraesCS1A02G065000.1    1020    LOC_Os05g28510.2        2607    80.682  176     29      4       722     893     2144    2318    1.36e-29        132
# TraesCS1A02G065000.1    1020    LOC_Os05g28510.4        2607    80.682  176     29      4       722     893     2144    2318    1.36e-29        132
# TraesCS1A02G065000.1    1020    LOC_Os05g28510.3        2607    80.682  176     29      4       722     893     2144    2318    1.36e-29        132
# TraesCS1A02G065000.1    1020    LOC_Os05g28510.1        2607    80.682  176     29      4       722     893     2144    2318    1.36e-29        132
# TraesCS1A02G065400.1    1986    LOC_Os01g28690.2        1854    78.995  438     88      4       1536    1971    1398    1833    8.95e-79        296
# TraesCS1A02G065400.1    1986    LOC_Os01g28690.1        1854    78.995  438     88      4       1536    1971    1398    1833    8.95e-79        296
# TraesCS1A02G065600.1    2610    LOC_Os05g28510.2        2607    85.402  2610    378     2       1       2610    1       2607    0.0     2706
# TraesCS1A02G065600.1    2610    LOC_Os05g28510.4        2607    85.402  2610    378     2       1       2610    1       2607    0.0     2706


def parse_blast(filepath, cutoff):
    wheat2rice = {}
    with open(filepath, 'r') as infh:
        for buf in infh:
            r = buf.replace('\n', '').split('\t')
            qid = r[0]         # query seq id length (wheat)
            qlen = int(r[1])   #    seq length
            sid = r[2]         # subject seq id (rice)
            slen = int(r[3])   #    seq length
            alen = int(r[5])   # alignment length
            
            if (1.0 * alen / qlen) > cutoff:
                if qid not in wheat2rice:
                    wheat2rice[qid] = []
                wheat2rice[qid].append(sid)
    return wheat2rice




def parse_go(filepath):
    rice2go = {}
    with open(filepath, 'r') as infh:
        for buf in infh:
            r = buf.replace('\n', '').split('\t')
            riceid = r[0]
            goid = r[1]
            if riceid not in rice2go:
                rice2go[riceid] = []
            rice2go[riceid].append(goid)
    return rice2go



def make_wheat_go(wheat2rice, rice2go, output_prefix):
    wheat2go = {}
    for wheatid, riceids in wheat2rice.items():
        for riceid in riceids:
            if riceid in rice2go:
                goterms = rice2go[riceid]
                if wheatid not in wheat2go:
                    wheat2go[wheatid] = []
                wheat2go[wheatid].extend(goterms)
    
    
    with open(output_prefix + '.cds.txt', 'w') as outfh:
        for wheatid in sorted(wheat2go.keys()):
            outfh.write('{0}\t{1}\n'.format(wheatid, ';'.join(sorted(list(set(wheat2go[wheatid]))))))
    
    gene2go = {}
    for wheatid, goterms in wheat2go.items():
        geneid, tsver = wheatid.split('.')
        if geneid not in gene2go:
            gene2go[geneid] = []
        gene2go[geneid].extend(goterms)
    
    with open(output_prefix + '.gene.txt', 'w') as outfh:
        for geneid in sorted(gene2go.keys()):
            outfh.write('{0}\t{1}\n'.format(geneid, ';'.join(sorted(list(set(gene2go[geneid]))))))
    
    

if __name__ == '__main__':
    
    blasttxt = sys.argv[1]
    ricego = sys.argv[2]
    cutoff = float(sys.argv[3])
    outprefix = sys.argv[4]
    
    wheat2rice = parse_blast(blasttxt, cutoff)
    rice2go = parse_go(ricego)
    make_wheat_go(wheat2rice, rice2go, outprefix)
    
    
    
