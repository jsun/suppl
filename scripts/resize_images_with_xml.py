import sys
import os
from pandabox import imgUtils

def resize_img(input_dpath, output_dpath, size):
    imgutils = imgUtils()
    imgutils.resize(input_dpath, (size, size), output_dpath, 'resized', n_jobs=-1)

if __name__ == '__main__':
    input_dpath = sys.argv[1]
    output_dpath = sys.argv[2]
    size = int(sys.argv[3])
    resize_img(input_dpath, output_dpath, size)

