import os
import sys
import argparse
from models import *



def predict(model_arch, model_path, class_labels, inference_dataset, mesh=None):
    
    dragonfly = DragonflyCls(model_arch=model_arch, model_path=model_path, class_labels=class_labels)
    probs = dragonfly.inference(inference_dataset)
    
    if mesh is not None:
        dragonflymesh = DragonflyMesh(mesh=mesh)
        mesh_output = dragonflymesh.inference(inference_dataset, k=2)
        probs = probs * mesh_output
    
    return probs
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')

    parser.add_argument('--class-label', required=True)
    parser.add_argument('--model-arch', required=True)
    parser.add_argument('--model-weight', required=True)
    parser.add_argument('--mesh', default=None)
    parser.add_argument('-i', '--inference-dataset', default=None)
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('--overwrite', action='store_true')
    
    args = parser.parse_args()
    
    probs = predict(args.model_arch, args.model_weight, args.class_label, args.inference_dataset, args.mesh)
    
    if args.output is None:
        print(probs)
    else:
        if args.overwrite and os.path.exists(args.output):
            probs.to_csv(args.output, header=False, index=True, sep='\t', mode='a')
        else:
            probs.to_csv(args.output, header=True, index=True, sep='\t')
    



