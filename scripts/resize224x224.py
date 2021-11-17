import sys
import os
import glob
import joblib
import cv2
import numpy as np



def resize(fpath):
    img_path = fpath['input']
    foutput = fpath['output']
    x = cv2.imread(img_path, cv2.IMREAD_COLOR)
    h, w, c = x.shape
    longest_edge = max(h, w)
    top = 0
    bottom = 0
    left = 0
    right = 0
    if h < longest_edge:
        diff_h = longest_edge - h
        top = diff_h // 2
        bottom = diff_h - top
    elif w < longest_edge:
        diff_w = longest_edge - w
        left = diff_w // 2
        right = diff_w - left
    elif h == longest_edge or w == longest_edge:
        pass
    else:
        print([h, w, c])
        raise ValueError('unexpected size.')
    x_resized = cv2.copyMakeBorder(x, top, bottom, left, right,
                           cv2.BORDER_CONSTANT, value=[0, 0, 0])
    x_resized = cv2.resize(x_resized, (224, 224))
    cv2.imwrite(foutput, x_resized)


def resize_images(input_dpath, output_dpath):
    
    image_files = []
    for input_fpath in glob.glob(os.path.join(input_dpath, '*', '*')):
        output_fpath = input_fpath.replace(input_dpath, output_dpath)
        image_files.append({'input': input_fpath, 'output': output_fpath})
        os.makedirs(os.path.dirname(output_fpath), exist_ok=True)
    
    
    joblib.Parallel(n_jobs=16, verbose=1)(joblib.delayed(resize)(image_file) for image_file in image_files)



if __name__ == '__main__':
    input_dpath = sys.argv[1]
    output_dpath = sys.argv[2]
    resize_images(input_dpath, output_dpath)



