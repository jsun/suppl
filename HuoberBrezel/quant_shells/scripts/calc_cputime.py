import sys

def get_time(x):

    y = x.replace('\n', '').split('\t')[1]

    y1, y2 = y.split('m')
    y2 = y2.replace('s', '')

    y1 = float(y1)
    y2 = float(y2)

    return y1 * 60 + y2



if __name__ == '__main__':

    infpath = sys.argv[1]
    nitems = int(sys.argv[2])  # HISAT, STAR : 1
                               # EAGLE-STAR: 3, EAGLE-COUNT: 6
    
    
    with open(infpath, 'r') as infh:
        
        i = 0
        s = 0
        for buf in infh:
            if buf[0:4] == 'real':
                buf = buf.replace('\n', '')
                i = i + 1
                s = s + get_time(buf)
                
                if i == nitems:
                    print(s)
                    i = 0
                    s = 0



