import os
import sys
import random
import numpy as np
import torch
import argparse
from model import *
import matplotlib.pyplot as plt


def train_cv(model, output, train_dataset, valid_dataset, batch_size=1024, epochs=100):
    
    if output is None:
        raise ValueError('argument `output` cannot be None.')
    
    
    # hyper-parameters of network architecture
    if model == 'L1':
        dropouts = [0, 0.2, 0.5, 0.6]
        n_hiddens = [[i] for i in [8, 10, 12, 14, 16, 18, 20, 22]]
        activate_funcs = ['relu', 'sigmoid']
        
    elif model == 'L2':
        dropouts = [0.5]
        n_hiddens = [[i, j] for i in [ 6, 10, 14, 18, 22, 24]
                            for j in [20, 24, 28, 32, 36, 40]]
        activate_funcs = ['relu']
        
    elif model == 'L3':
        dropouts = [0.5]
        n_hiddens = [[i, j, k] for i in [ 6, 10, 14, 18, 22, 24]
                               for j in [20, 24, 28, 32, 36, 40]
                               for k in [20, 24, 28, 32, 36, 40]]
        activate_funcs = ['relu']
    
    
    eval_stats = {'activate_fun': [], 'dropout': [], 'train_mse': [], 'valid_mse': []}
    for i in range(len(n_hiddens[0])):
        eval_stats['n_hidden_' + str(i)] = []
    
    
    for activate_func in activate_funcs:
        for n_hidden in n_hiddens:
            for dropout in dropouts:
                seed = int(dropout * 100 + len(n_hidden) * n_hidden[0] + n_hidden[int((len(n_hidden) - 1)/2)] * n_hidden[-1])
                random.seed(seed)
                np.random.seed(seed)
                torch.manual_seed(seed)
                
                output_ = '{}_{}_{}_{}'.format(output, activate_func, '-'.join([str(i) for i in n_hidden]), dropout)
                loss_hist = train(output_, train_dataset, valid_dataset, batch_size, epochs,
                                  n_hidden, dropout, activate_func)
                
                n_tries = len(loss_hist['train'])
                
                eval_stats['activate_fun'].extend([activate_func] * n_tries)
                eval_stats['dropout'].extend([dropout] * n_tries)
                eval_stats['train_mse'].extend(loss_hist['train'])
                eval_stats['valid_mse'].extend(loss_hist['valid'])
                for i in range(len(n_hidden)):
                    eval_stats['n_hidden_' + str(i)] = n_hidden[i]
                
                pd.DataFrame(eval_stats).to_csv(output + '_cv_stats.tsv', index=False, header=True, sep='\t')




    
def train(output, train_dataset, valid_dataset, batch_size, epochs, n_hidden=16, dropout=0.5, activate_func='relu'):
    
    if isinstance(n_hidden, int):
            n_hidden = [n_hidden]
    
    if output is None:
        raise ValueError('argument `output` cannot be None.')
    
    
    # validate 10 times (dropout the initialization effects)
    min_loss = {'train': [], 'valid': []}
    
    for i in range(10):
        jppnet = Jppnet(n_hidden=n_hidden, dropout=dropout, activate_func=activate_func)
        train_history = jppnet.train(train_dataset, valid_dataset, batch_size, epochs)

        min_loss['train'].append(train_history.loc[:, 'train'].iloc[-1])#.min())
        min_loss['valid'].append(train_history.loc[:, 'valid'].iloc[-1])#.min())
    
    print('-------------')
    print('n_hidden: {}  dropout: {}  activate_func: {}'.format('-'.join([str(i) for i in n_hidden]), dropout, activate_func))
    print('train loss: {};    valid loss: {}'.format(np.mean(min_loss['train']),
                                                     np.mean(min_loss['valid'])))
    
    #jppnet.save(output)
    #pd.DataFrame(min_loss).to_csv(output + '_stats.tsv', index=False, header=True, sep='\t')
    
    return min_loss
    
    



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')
     
    parser.add_argument('--model', default='L1')
    parser.add_argument('--output', default='./weights/test.pth')
    parser.add_argument('--train-dataset', default='../formatted_data/kyuuri_honpo_percent.train.tsv')
    parser.add_argument('--valid-dataset', default='../formatted_data/kyuuri_honpo_percent.valid.tsv')
    parser.add_argument('--epochs', default=100, type=int)
    parser.add_argument('--batch-size', default=1024, type=int)
    parser.add_argument('--train-mode', default='train')
    
    args = parser.parse_args()
    
    if args.train_mode == 'cv':
        train_cv(args.model, args.output, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs)
        
    elif args.train_mode == 'train':
        train(args.output, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs)
    



'''
python train.py -t ../formatted_data/kyuuri_honpo_percent.all.tsv \
                -v ../formatted_data/kyuuri_honpo_percent.all.tsv


'''


