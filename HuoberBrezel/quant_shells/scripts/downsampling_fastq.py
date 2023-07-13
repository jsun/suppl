import os
import sys
import gzip
import random


def downsampling_fasq(input_fpath, output_fpath, downsample_rate):
    
    random.seed(int(hash(os.path.basename(input_fpath))))
    
    record = []
    i = 0

    with gzip.open(input_fpath, 'rb') as infh, gzip.open(output_fpath, 'wb') as outfh:

        for buf in infh:
            i = i + 1

            record.append(buf)

            if i == 4:
                if random.random() <= downsample_rate:
                    for r in record:
                        outfh.write(r)

                i = 0
                record = []




if __name__ == '__main__':

    input_fpath = sys.argv[1]
    output_fpath = sys.argv[2]
    downsample_rate = float(sys.argv[3])

    downsampling_fasq(input_fpath, output_fpath, downsample_rate)



