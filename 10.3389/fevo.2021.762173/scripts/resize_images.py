import sys
import cv2
import numpy as np



def resize(img_path, foutput = None, w = 1024):
    resized = 0
    
    img = cv2.imread(img_path, cv2.IMREAD_COLOR)
    img_h, img_w = img.shape[:2]
    if img_w > w:
        r = w / img_w

        h = int(img_h * r)

        img_resized = cv2.resize(img, (w, h))
        if foutput.rfind('.jpg') > 0:
            cv2.imwrite(foutput, img_resized, [int(cv2.IMWRITE_JPEG_QUALITY), 100])
        else:
            cv2.imwrite(foutput, img_resized)
        
        resized = 1
    
    return resized
    

if __name__ == '__main__':

    finput = sys.argv[1]
    foutput = sys.argv[2]

    resized_signal = resize(finput, foutput)
    exit(resized_signal)



