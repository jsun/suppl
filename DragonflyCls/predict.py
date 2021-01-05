import os
import sys
import argparse
from models import *



PROJECT_PATH = os.path.dirname(__file__)
IMAGE_FOLDER = os.path.join(PROJECT_PATH, 'testimagedir/Sympetrum_darwinianum___4.cropped.png')
WEIGHTS_PATH = os.path.join(PROJECT_PATH, 'weights/dragonfly-cls-model.pth')
CLASS_PATH = os.path.join(PROJECT_PATH, 'data/class_labels.txt')


def predict(model_arch, model_path, class_labels, inference_dataset):
    
    dragonfly = DragonflyCls(model_arch=model_arch, model_path=model_path, class_labels=class_labels)
    ret = dragonfly.inference(inference_dataset)
    
    return ret
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')

    parser.add_argument('--class-label', required=True)
    parser.add_argument('--model-arch', required=True)
    parser.add_argument('--model-weight', required=True)
    parser.add_argument('-i', '--inference-dataset', default=IMAGE_FOLDER)
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()
    
    ret = predict(args.model_arch, args.model_weight, args.class_label, args.inference_dataset)
    if args.output is None:
        print(ret)
    else:
        ret.to_csv(args.output, header=True, index=True, sep='\t')
    



