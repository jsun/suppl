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
    dropouts = [0, 0.5]
    n_hiddens_1 = n_hiddens_2 = [8, 10, 12, 14, 16, 18, 20, 22]
    activate_funcs = ['relu', 'sigmoid']
    
    n_zeros = len(dropouts) * len(n_hiddens_1) * len(n_hiddens_2) * len(activate_funcs)
    eval_stats = {'activate_fun': [], 'n_hidden_1': [], 'n_hidden_2': [],
                  'dropout': [], 'train_mse': [], 'valid_mse': []}
    
    
    for activate_func in activate_funcs:
        for n_hidden_1 in n_hiddens_1:
            for n_hidden_2 in n_hiddens_2:
                for dropout in dropouts:
    
                    random.seed(int(dropout * 100 + n_hidden_1 + n_hidden_2 * 10))
                    np.random.seed(int(dropout * 100 + n_hidden_1 + n_hidden_2 * 10))
                    torch.manual_seed(int(dropout * 100 + n_hidden_1 + n_hidden_2 * 10))
                    
                    output_ = '{}_{}_{}_{}_{}'.format(output, activate_func,
                                 n_hidden_1, n_hidden_2, dropout)
                    loss_hist = train(model, output_, train_dataset, valid_dataset, batch_size, epochs,
                                  n_hidden_1, n_hidden_2, dropout, activate_func)
                
                    n_tries = len(loss_hist['train'])
              
                    eval_stats['activate_fun'].extend([activate_func] * n_tries)
                    eval_stats['n_hidden_1'].extend([n_hidden_1] * n_tries)
                    eval_stats['n_hidden_2'].extend([n_hidden_2] * n_tries)
                    eval_stats['dropout'].extend([dropout] * n_tries)
                    eval_stats['train_mse'].extend(loss_hist['train'])
                    eval_stats['valid_mse'].extend(loss_hist['valid'])
    
                    # save every tries
                    pd.DataFrame(eval_stats).to_csv(output + '_cv_stats.tsv', index=False, header=True, sep='\t')




    
def train(model, output, train_dataset, valid_dataset, batch_size, epochs,
            n_hidden_1=16, n_hidden_2=16, dropout=0.5, activate_func='relu'):
    
    if output is None:
        raise ValueError('argument `output` cannot be None.')
    
    
    # validate 10 times (dropout the initialization effects)
    min_loss = {'train': [], 'valid': []}
    
    for i in range(10):
        jppnet = Jppnet(model_arch=model, n_hidden_1=n_hidden_1, n_hidden_2=n_hidden_2, 
                        dropout=dropout, activate_func=activate_func)
        train_history = jppnet.train(train_dataset, valid_dataset, batch_size, epochs)

        min_loss['train'].append(train_history.loc[:, 'train'].min())
        min_loss['valid'].append(train_history.loc[:, 'valid'].min())
    
    print('-------------')
    print('n_hidden_1: {}  n_hidden_2: {}    dropout: {}  activate_func: {}'.format(n_hidden_1, n_hidden_2, dropout, activate_func))
    print('train loss: {};    valid loss: {}'.format(np.mean(min_loss['train']),
                                                     np.mean(min_loss['valid'])))
    
    jppnet.save(output)
    pd.DataFrame(min_loss).to_csv(output + '_stats.tsv', index=False, header=True, sep='\t')
    
    return min_loss
    
    



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')
     
    parser.add_argument('--model', default='2L')
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
        train(args.model, args.output, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs)
    



'''
python train.py -t ../formatted_data/kyuuri_honpo_percent.all.tsv \
                -v ../formatted_data/kyuuri_honpo_percent.all.tsv


'''


