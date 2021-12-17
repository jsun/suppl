import os
import sys


in_fname = sys.argv[1]

out_fname = in_fname.replace(' ', 'sp')
out_fname = out_fname.replace('(', 'lb')
out_fname = out_fname.replace(')', 'rb')

os.rename(in_fname, out_fname)

