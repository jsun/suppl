import os
import sys
import random
import json
import numpy as np
import torch
import argparse
from model import *
import matplotlib.pyplot as plt
import sklearn.metrics


TIRAMISU_DEBUG = False


def get_params(params_fpath):
    
    if params_fpath is None:
        params = {
            'dropout': 0.5,
            #'n_hidden': [10, 32, 32],
            'n_hidden': [21, 40],
            'activate_func': 'relu'
        }
    
    else:
        with open(params_fpath, 'r') as jsonfh:
            params = json.load(jsonfh)
        
        n_hidden = []
        if 'n_hidden_0' in params:
            n_hidden.append(params['n_hidden_0'])
            del params['n_hidden_0']
        if 'n_hidden_1' in params:
            n_hidden.append(params['n_hidden_1'])
            del params['n_hidden_1']
        if 'n_hidden_2' in params:
            n_hidden.append(params['n_hidden_2'])
            del params['n_hidden_2']
    
        params['n_hidden'] = n_hidden
    
    return params
    



def train_cv(model, weight, train_dataset, valid_dataset, batch_size=1024, epochs=50, n_tries=1):
    
    if weight is None:
        raise ValueError('argument `weight` cannot be None.')
    
    # hyper-parameters of network architecture
    if model == 'L1':
        dropouts = [0, 0.5]
        n_hiddens = [[i] for i in [4, 8, 12, 16, 20, 24, 28]]
        activate_funcs = ['relu']
        if TIRAMISU_DEBUG:
            dropouts = [0, 0.5]
            n_hiddens = [[i] for i in [8, 22]]
            activate_funcs = ['relu']
        
    elif model == 'L2':
        dropouts = [0.5]
        n_hiddens = [[i, j] for i in [8, 16, 24, 32, 40]
                            for j in [8, 16, 24, 32, 40]]
        activate_funcs = ['relu']
        if TIRAMISU_DEBUG:
            dropouts = [0, 0.5]
            n_hiddens = [[i, j] for i in [ 8, 22]
                                for j in [18, 24]]
            activate_funcs = ['relu']
        
    elif model == 'L3':
        dropouts = [0.5]
        n_hiddens = [[i, j, k] for i in [ 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
                               for j in [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]
                               for k in [20, 24, 28, 32, 36, 40]]
        activate_funcs = ['relu']
    
    
    eval_stats = {'activate_func': [], 'dropout': [], 'train_mse': [], 'valid_mse': []}
    for i in range(len(n_hiddens[0])):
        eval_stats['n_hidden_' + str(i)] = []
    
    for activate_func in activate_funcs:
        for n_hidden in n_hiddens:
            for dropout in dropouts:
                seed = int(dropout * 100 + len(n_hidden) * n_hidden[0] + n_hidden[int((len(n_hidden) - 1)/2)] * n_hidden[-1])
                random.seed(seed)
                np.random.seed(seed)
                torch.manual_seed(seed)
                
                weight_ = '{}_{}_{}_{}'.format(weight, activate_func, '-'.join([str(i) for i in n_hidden]), dropout)
                param_ = {'dropout': dropout, 'n_hidden': n_hidden, 'activate_func': activate_func}
                print('--------')
                print(param_)
                loss_hist = train(weight_, train_dataset, valid_dataset, batch_size, epochs, param_, False, n_tries)
                
                eval_stats['activate_func'].extend([activate_func] * n_tries)
                eval_stats['dropout'].extend([dropout] * n_tries)
                eval_stats['train_mse'].extend(loss_hist['train'])
                eval_stats['valid_mse'].extend(loss_hist['valid'])
                for i in range(len(n_hidden)):
                    eval_stats['n_hidden_' + str(i)].extend([n_hidden[i]] * n_tries)
                
                pd.DataFrame(eval_stats).to_csv(weight + '_cv_stats.tsv', index=False, header=True, sep='\t')




    
def train(weight, train_dataset, valid_dataset, batch_size, epochs, params = None, save_weight=True, n_try=10):
    if TIRAMISU_DEBUG:
        epochs = 2
    
    if params is None:
        params = get_params()
    
    if isinstance(params['n_hidden'], int):
            params['n_hidden'] = [params['n_hidden']]
    
    if weight is None:
        raise ValueError('argument `weight` cannot be None.')
    
    
    best_loss_ = 1e10
    best_model_ = None
    min_loss = {'train': [], 'valid': []}
    train_history_ = None
    
    for i in range(n_try):
        jppnet = Jppnet(n_hidden=params['n_hidden'], dropout=params['dropout'], activate_func=params['activate_func'])
        train_history = jppnet.train(train_dataset, valid_dataset, batch_size, epochs)
        
        min_loss['train'].append(train_history.loc[:, 'train'].min())
        min_loss['valid'].append(train_history.loc[:, 'valid'].min())
        
        if train_history.loc[:, 'valid'].min() < best_loss_:
            best_loss_  = train_history.loc[:, 'valid'].min()
            best_model_ = jppnet
            train_history_ = train_history
    
    print('train loss: {};    valid loss: {}'.format(np.mean(min_loss['train']),
                                                     np.mean(min_loss['valid'])))
    
    if save_weight:
        jppnet.save(weight)
        pd.DataFrame(min_loss).to_csv(weight + '_stats.tsv', index=False, header=True, sep='\t')
        pd.DataFrame(train_history_).to_csv(weight + '_trainhistory.tsv', index=False, header=True, sep='\t')
    
    return min_loss




def validate(weight, params, valid_dataset, output):
    
    jppnet = Jppnet(n_hidden=params['n_hidden'], dropout=params['dropout'], activate_func=params['activate_func'])
    jppnet.load_weights(weight)
    predicted_df = jppnet.validate(valid_dataset)
    predicted_df.to_csv(output + '.tsv', index=False, header=True, sep='\t')
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.scatter(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])
    ax.set_xlabel('label')
    ax.set_ylabel('predicted')
    fig.savefig(output + '.png')
    
    print('--- validation results ---')
    print('R2 :   {}'.format(sklearn.metrics.r2_score(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('MSE:   {}'.format(sklearn.metrics.mean_squared_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('RMSE:  {}'.format(np.sqrt(sklearn.metrics.mean_squared_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted']))))
    print('MAE:   {}'.format(sklearn.metrics.mean_absolute_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('--------------------------')
    
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')
     
    parser.add_argument('--model', default='L1')
    parser.add_argument('--weight', default='./weights/test.pth')
    parser.add_argument('--output', default='./weights/test.txt')
    parser.add_argument('--train-dataset', default='../formatted_data/test/train_std.tsv')
    parser.add_argument('--valid-dataset', default='../formatted_data/test/valid_std.tsv')
    parser.add_argument('--epochs', default=50, type=int)
    parser.add_argument('--batch-size', default=1024, type=int)
    parser.add_argument('--mode', default='train')
    parser.add_argument('--params', default=None)
    
    args = parser.parse_args()
    
    if args.mode == 'cv':
        train_cv(args.model, args.weight, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs)
        
    elif args.mode == 'train':
        params = get_params(args.params)
        train(args.weight, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs, params, True, 3)
    
    elif args.mode == 'valid':
        params = get_params(args.params)
        validate(args.weight, params, args.valid_dataset, args.output)
    
    elif args.mode == 'inference':
        pass    
    
    
    


