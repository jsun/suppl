with open('AT_GO_all_subset.txt', 'r') as fh:
    for buf in fh:
        buf = buf.replace('\n', '')
        buf_record = buf.split('\t')
        go_record = list(set(buf_record[1].split(',')))
        
        #if (10 <= len(go_record)) and (len(go_record) <= 500):
        print(buf_record[0] + '\t' + ','.join(go_record))

            
