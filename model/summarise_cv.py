import os
import sys
import glob
import json
import pandas as pd



def save_best_params(model_arch, cv_dpath, output_fpath):
    
    stats_df = None
    
    for fpath in glob.glob(os.path.join(cv_dpath, 'weight_' + model_arch + '.cv_*.tsv')):
        stats_df_ = pd.read_csv(fpath, sep='\t', header=0)
        if stats_df is None:
            stats_df = stats_df_
        else:
            stats_df = pd.concat([stats_df, stats_df_])
    
    if model_arch == 'L1':
        stats_df = stats_df.groupby(['activate_func', 'dropout', 'n_hidden_0']).mean().reset_index()
    elif model_arch == 'L2':
        stats_df = stats_df.groupby(['activate_func', 'dropout', 'n_hidden_0', 'n_hidden_1']).mean().reset_index()
    
    min_mse_idx = int(stats_df[['valid_mse']].idxmin())
    stats_df.iloc[min_mse_idx, :].drop(labels = ['train_mse', 'valid_mse']).to_json(output_fpath)



if __name__ == '__main__':
    
    args =sys.argv
    model_arch = args[1]
    cv_dpath = args[2]
    output_fpath = args[3]
    
    save_best_params(model_arch, cv_dpath, output_fpath)


