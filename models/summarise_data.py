import os
import sys
import glob
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
sns.set_style('whitegrid')
sns.set_palette('Greys_r')



def plot_hist(dpath):
    x = pd.read_csv(os.path.join(dpath, 'data.tsv'), sep='\t', header=0)
    x_incidence = x.loc[:, 'incidence']
    x_incidence = x_incidence.dropna()
    
    fpath = os.path.join(dpath, 'data_distr_hist.png')

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.hist(x_incidence, bins = 20)
    fig.savefig(fpath)
    


if __name__ == '__main__':
    
    args =sys.argv
    dpath = args[1]
    
    plot_hist(dpath)
    
    
    
    

