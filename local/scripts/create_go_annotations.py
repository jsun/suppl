import os
import sys
import re
import argparse

#
# ptyhon $this --obo go.obo --goslim ATH_GO_GOSLIM.txt --ann chi_m25.txt
#

def get_parental_goterms(goterm, go, del_targets):
    if goterm in go:
        parental_goterms = go[goterm]
        del_targets.extend(parental_goterms)
        
        parental_parental_goterms = []
        for i in range(len(parental_goterms)):
            _parental_parental_goterms = get_parental_goterms(parental_goterms[i], go, del_targets)
            
        


def extend_list(l):
    j = []
    for i in l:
        if type(i) is list:
            k = extend_list(i)
            j.extend(k)
        else:
            j.append(i)
    return j


def calc_go_levels(goterm, go, go_level):
    if goterm not in go:
        return go_level
    else:
        parental_goterms = go[goterm] 
        go_levels = [go_level + 1] * len(parental_goterms)
        
        for i in range(len(parental_goterms)):
            go_levels[i] = calc_go_levels(parental_goterms[i], go, go_levels[i])
        
        return go_levels
        
        

def proc_main(obofile, goslimfile, annfile, go2carhrfile):
    
    ## create ID mapping object for exchanging ID from TAIR to CARHR.
    tair2carhr = {}
    with open(annfile, 'r') as anfh:
        for buf in anfh:
            buf_records = buf.split('\t')
            tair2carhr[buf_records[1]] = buf_records[0]
    
    ## replace TAIR GOSlim data with CARHR ID.
    carhr2go = {}
    go2carhr = {}
    with open(goslimfile, 'r') as gofh:
        for buf in gofh:
            buf_records = buf.split('\t')
            if buf_records[0] in tair2carhr:
                ## 'CARHR282590': ['GO:0003674', 'GO:0008150']
                if tair2carhr[buf_records[0]] not in carhr2go:
                    carhr2go[tair2carhr[buf_records[0]]] = []
                carhr2go[tair2carhr[buf_records[0]]].append(buf_records[5])
                ## 'GO:0000150': ['CARHR121570', 'CARHR140490']
                if buf_records[5] not in go2carhr:
                    go2carhr[buf_records[5]] = []
                go2carhr[buf_records[5]].append(tair2carhr[buf_records[0]])
    
    for goterm in go2carhr.keys():
        go2carhr[goterm] = list(set(go2carhr[goterm]))
    
   
    ## create GO object for calculating the levels
    go = {}
    with open(obofile, 'r') as obofh:
        for buf in obofh:
            if buf[0:6] == '[Term]':
                go_term = ''
                isBP = False
            if buf[0:4] == 'id: ':
                go_term = buf[4:14]
            if buf[0:11] == 'namespace: ':
                if 'biological_process' in buf:
                    isBP = True
            if buf[0:6] == 'is_a: ':
                if isBP:
                    if go_term not in go:
                        go[go_term] = []
                    go[go_term].append(buf[6:16])
    
    
    ## delete GO terms higher at 1-3 levels.
    go_under4 = {}
    for goterm in go.keys():
        goterm_lower_level = max(extend_list(calc_go_levels(goterm, go, go_level = 1)))
        if goterm_lower_level > 3:
            go_under4[goterm] = go[goterm]
    go = go_under4
    

    # remove the parental terms.
    del_targets = []
    for goterm in go.keys():
        get_parental_goterms(goterm, go, del_targets)
    del_targets = list(set(del_targets))
    
    ## delete GO terms (consists of <10 or >500 genes) and parental temrs.
    #go2carhr_filtered = {}
    #carhr2go_filtered = {}
    #for goterm in go2carhr.keys():
    #    if goterm in go:
    #        if 10 < len(go2carhr[goterm]) and len(go2carhr[goterm]) < 500:
    #            if goterm not in del_targets:
    #                go2carhr_filtered[goterm] = go2carhr[goterm]
    #                for gene in go2carhr_filtered[goterm]:
    #                    if gene not in carhr2go_filtered:
    #                        carhr2go_filtered[gene] = []
    #                    carhr2go_filtered[gene].append(goterm)
    
    
    ## delete GO terms (consists of <10 or >500 temrs)
    go2carhr_filtered = {}
    carhr2go_filtered = {}
    for goterm in go2carhr.keys():
        if goterm in go:
            if 10 < len(go2carhr[goterm]) and len(go2carhr[goterm]) < 500:
                    go2carhr_filtered[goterm] = go2carhr[goterm];
                
    
    
    with open(go2carhrfile, 'w') as go2carhrfh:
        for goterm in go2carhr_filtered.keys():
            go2carhrfh.write(goterm + '\t')
            go2carhrfh.write(';'.join(list(set(go2carhr_filtered[goterm]))) + '\n')
    
    
    #with open(carhr2gofile, 'w') as carhr2gofh:
    #    for gene in carhr2go_filtered.keys():
    #        carhr2gofh.write(gene + '\t')
    #        carhr2gofh.write(';'.join(list(set(carhr2go_filtered[gene]))) + '\n')
    




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'GO data generation.')
    parser.add_argument('-o', '--obo', required = True)
    parser.add_argument('-g', '--goslim', required = True)
    parser.add_argument('-a', '--ann', required = True)
    parser.add_argument('-j', '--go2carhr', required = True)
    args = parser.parse_args()
    
    proc_main(args.obo, args.goslim, args.ann, args.go2carhr)
    
