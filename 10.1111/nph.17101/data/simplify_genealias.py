

tair2name = {}
tair2fullname = {}

with open('gene_aliases_20140331.txt', 'r') as infh:
    for buf in infh:
        if buf[0:10] == 'locus_name':
            continue
        tair_id, tair_name, tair_fullname = buf.replace('\n', '').replace('#', '').replace('"', '').replace('\'', '').split('\t')
            
        if tair_id not in tair2name:
            tair2name[tair_id] = set()
            tair2fullname[tair_id] = set()
        
        if tair_name != '':
            tair2name[tair_id].add(tair_name)
        if tair_fullname != '':
            tair2fullname[tair_id].add(tair_fullname)


for tair_id in tair2name.keys():
    print(tair_id + '\t' + ' '.join(tair2name[tair_id]) + '\t' + '; '.join(tair2fullname[tair_id]))
    
