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
sns.set_palette('Set1')


def calc_rmse(y_obs, y_est):
    y_obs = np.array(y_obs)
    y_est = np.array(y_est)
    y_est[y_est < 0] = 0
    y_diff = y_obs - y_est
    rmse = np.sqrt(np.sum(y_diff ** 2) / len(y_obs))
    
    return rmse



def plot_hist(rmse, fpath):
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.hist([rmse['norm_datasets'], rmse['null_datasets']], bins = 20, stacked=False,
            label=['normal', 'random'])
    ax.legend()
    fig.savefig(fpath)
    
    


def calc_validscore(dpath):
    
    rmse = {'norm_datasets': [], 'null_datasets': []}
    
    for ds in rmse.keys():
        for dsid in range(100):
            dsid = '{:03d}'.format(dsid + 1)
            fpath = os.path.join(dpath, ds, dsid, 'modelvalid_L2.tsv')
            
            if os.path.exists(fpath):
                dat = pd.read_csv(fpath, sep='\t', header=0)
                rmse[ds].append(calc_rmse(dat.iloc[:, 0], dat.iloc[:, 1]))
    
    pd.DataFrame(rmse).to_csv('testfile.csv') 
              
    fpath = os.path.join(dpath, 'rmse_hist.png')
    plot_hist(rmse, fpath)
    
    norm_mu = np.mean(rmse['norm_datasets'])
    norm_sd = np.std(rmse['norm_datasets'])
    null_mu = np.mean(rmse['null_datasets'])
    null_sd = np.std(rmse['null_datasets'])
    z_score = (norm_mu - null_mu) / null_sd
    
    stats = [os.path.basename(dpath), str(norm_mu), str(norm_sd), str(null_mu), str(null_sd), str(z_score)]
    stats = '\t'.join(stats)
    print(stats)


if __name__ == '__main__':
    
    args =sys.argv
    dpath = args[1]
    
    calc_validscore(dpath)
    
    
    
    

