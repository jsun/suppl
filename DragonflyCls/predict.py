import os
import sys
import argparse
from models import *



PROJECT_PATH = os.path.dirname(__file__)
IMAGE_FOLDER = os.path.join(PROJECT_PATH, 'testimagedir/Sympetrum_darwinianum___4.cropped.png')
WEIGHTS_PATH = os.path.join(PROJECT_PATH, 'weights/dragonfly-cls-model.pth')
CLASS_PATH = os.path.join(PROJECT_PATH, 'data/class_labels.txt')





def predict(inference_dataset, model=WEIGHTS_PATH, class_label=CLASS_PATH):
    
    dragonfly = DragonflyCls(class_labels=class_label, model=model)
    ret = dragonfly.inference(inference_dataset)
    
    return ret
    



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')

    parser.add_argument('-m', '--model', default=WEIGHTS_PATH)
    parser.add_argument('-c', '--class-label', default=CLASS_PATH)
    parser.add_argument('-i', '--inference-dataset', default=IMAGE_FOLDER)
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()
    
    ret = predict(args.inference_dataset, args.model, args.class_label)
    if args.output is None:
        print(ret)
    else:
        ret.to_csv(args.output, header=True, index=True, sep='\t')
    



