
h = {}

with open('ATH_GO_GOSLIM.20180605.txt', 'r') as infh:
    for buf in infh:
        r = buf.split('\t')
        #if (r[7] == 'P'):
        if r[0] not in h:
            h[r[0]] = set()
        h[r[0]].add(r[5])
            

for k, v in h.items():
    if len(v) > 0 and len(v) < 500:
        print(k + '\t' + ','.join(v))


