import os
import sys
import argparse
import cv2
from models import *



def predict(model_arch, model_path, class_labels, inference_dataset, mesh=None, d=100):
    
    dragonfly = DragonflyCls(model_arch=model_arch, model_path=model_path, class_labels=class_labels, device='cpu')
    probs = dragonfly.inference(inference_dataset)
    
    if mesh is not None:
        dragonflymesh = DragonflyMesh(mesh=mesh)
        mesh_output = dragonflymesh.inference(inference_dataset, d=d)
        probs = probs * mesh_output
    
    return probs
    

def predict_gradcam(model_arch, model_path, class_labels, inference_dataset, output_path):
    
    input_images = []
    if os.path.isfile(inference_dataset):
        input_images.append(inference_dataset)
    else:
        for fpath in os.listdir(inference_dataset):
            if os.path.splitext(fpath)[1].lower() in ['.jpg', '.jpeg', '.png']:
                input_images.append(os.path.join(inference_dataset, fpath))
    
    dragonfly = DragonflyCls(model_arch=model_arch, model_path=model_path, class_labels=class_labels, device='cpu')
    
    for input_image in input_images:
        output_fpath = os.path.join(output_path, os.path.splitext(os.path.basename(input_image))[0] + '.gradcam.png')
        gradcam_img = dragonfly.gradcam(input_image)
        cv2.imwrite(output_fpath, gradcam_img)


    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')

    parser.add_argument('--class-label', required=True)
    parser.add_argument('--model-arch', required=True)
    parser.add_argument('--model-weight', required=True)
    parser.add_argument('--mesh', default=None)
    parser.add_argument('-d', default=None, type=int)
    parser.add_argument('-i', '--inference-dataset', default=None)
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-g', '--gradcam-output', default=None)
    parser.add_argument('--overwrite', action='store_true')
    
    args = parser.parse_args()
    
    if args.gradcam_output is not None:
        predict_gradcam(args.model_arch, args.model_weight, args.class_label,
                        args.inference_dataset, args.gradcam_output)
    
    else:
        probs = predict(args.model_arch, args.model_weight, args.class_label,
                        args.inference_dataset, args.mesh, args.d)
    
        if args.output is None:
            print(probs)
        else:
            if args.overwrite and os.path.exists(args.output):
                probs.to_csv(args.output, header=False, index=True, sep='\t', mode='a')
            else:
                probs.to_csv(args.output, header=True, index=True, sep='\t')


