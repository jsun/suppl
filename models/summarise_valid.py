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





def calc_rmse(y_obs, y_est):
    y_obs = np.array(y_obs)
    y_est = np.array(y_est)
    y_est[y_est < 0] = 0
    y_diff = y_obs - y_est
    rmse = np.sqrt(np.sum(y_diff ** 2) / len(y_obs))
    return rmse


def calc_rmse_from_df(df):
    return calc_rmse(df.loc[:, 'true'].values, df.loc[:, 'predicted'].values)


def summarise_classic_cvresults(dpath):
    sum_table = []
    
    for fpath in sorted(glob.glob(os.path.join(dpath, '*.tsv'))):
        fname = os.path.basename(fpath)
        dat = os.path.splitext(fname)[0].split('____')
        
        cv_results = pd.read_csv(fpath, sep='\t', header=0)
        rmse = []
        for k in cv_results.loc[:, 'cv'].unique():
            rmse.append(calc_rmse_from_df(cv_results.loc[cv_results.loc[:, 'cv'] == k, :]))
        
        dat.extend([np.mean(rmse), np.var(rmse, ddof=1)])
        dat.extend(rmse)
        
        sum_table.append(dat)
    
    sum_table = pd.DataFrame(sum_table,
                    columns=['model', 'feature_type', 'data_type', 'disease',
                             'cv_mean', 'cv_var',
                             'cv_1', 'cv_2', 'cv_3', 'cv_4', 'cv_5',
                             'cv_6', 'cv_7', 'cv_8', 'cv_9', 'cv_10'])
    return sum_table
    
    
def summarise_dnn_cvresults(dpath):
    sum_table = []
    
    for ddpath in sorted(glob.glob(os.path.join(dpath, '**'))):
        ddname = os.path.basename(ddpath)
        for fpath in sorted(glob.glob(os.path.join(ddpath, '*_cvresults.tsv'))):
            fname = os.path.basename(fpath)
            dat = fname.replace('_cvresults.tsv', '').split('____')
            
            cv_results = pd.read_csv(fpath, sep='\t', header=0)
            rmse = []
            for k in cv_results.loc[:, 'cv'].unique():
                rmse.append(calc_rmse_from_df(cv_results.loc[cv_results.loc[:, 'cv'] == k, :]))
        
            dat.extend([np.mean(rmse), np.var(rmse, ddof=1)])
            dat.extend(rmse)
            
            sum_table.append(dat)

    sum_table = pd.DataFrame(sum_table,
                    columns=['model', 'feature_type', 'data_type', 'disease',
                             'cv_mean', 'cv_var',
                             'cv_1', 'cv_2', 'cv_3', 'cv_4', 'cv_5',
                             'cv_6', 'cv_7', 'cv_8', 'cv_9', 'cv_10'])
    return sum_table
    
    

def main(dpath):
    sum_table_rc = summarise_classic_cvresults(os.path.join(dpath, 'randomize', 'classic_models'))
    sum_table_rd = summarise_dnn_cvresults(os.path.join(dpath, 'randomize', 'dnn_models'))
    sum_table_sc = summarise_classic_cvresults(os.path.join(dpath, 'shuffle', 'classic_models'))
    sum_table_sd = summarise_dnn_cvresults(os.path.join(dpath, 'shuffle', 'dnn_models'))
    
    sum_table = pd.concat([sum_table_rc, sum_table_rd, sum_table_sc, sum_table_sd])
    sum_table = sum_table.sort_values(['model', 'feature_type', 'disease', 'data_type']).reset_index(drop=True)
    
    sum_table.to_csv(os.path.join(dpath, 'summary.tsv'), sep='\t', header=True, index=False)
    

if __name__ == '__main__':
    
    args =sys.argv
    dpath = args[1]
    main(dpath)
    
    
    
    

