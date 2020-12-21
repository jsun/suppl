#python calc.go.num.py  > go.num.txt

go2tair = {}

with open('ATH_GO_GOSLIM.txt', 'r') as infh:
    for buf in infh:
        buf = buf.replace('\n', '')
        r = buf.split('\t')
        for goid in r[1].split(','):
            if goid not in go2tair:
                go2tair[goid] = set()
            go2tair[goid].add(r[0])



for goid, tairs in go2tair.items():
    print(goid + '\t' + str(len(tairs)))


