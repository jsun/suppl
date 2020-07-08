import os
import sys
import argparse
from model import *


def test(test_dataset, trained_model, output=None):
    
    
    dragonfly = miniDragonfly()
    dragonfly.load_weights(trained_model)
    pred = dragonfly.inference(test_dataset)
    
    if output is None:
        print(pred)
    else:
        pred.to_csv(output, index=False, sep='\t')
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')
   
    parser.add_argument('-m', '--model', default='./weights/jppnet-nn-model.pth')
    parser.add_argument('-t', '--test-dataset', default='../formatted_data/kyuuri_honpo_percent.test.tsv')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-b', '--batch-size', type=int, default=32)
    args = parser.parse_args()
    
    test(args.test_dataset, args.model, args.output)
    
    


'''
python test.py -t ../formatted_data/summer_kyuuri_honpo_percent.test.tsv \
               -o ../test.result.percente.txt 

'''
