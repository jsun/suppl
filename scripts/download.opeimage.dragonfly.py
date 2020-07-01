import os
import sys
import time
import urllib.request
import cv2


xml_header = '''\
<annotation>
        <folder>dragonfly</folder>
        <filename>{}</filename>
        <path>{}</path>
        <source>
                <database>dragonfly</database>
        </source>
        <size>
                <width>{}</width>
                <height>{}</height>
                <depth>3</depth>
        </size>
        <segmented>0</segmented>
'''

xml_record = '''\
        <object>
                <name>dragonfly</name>
                <pose>Unspecified</pose>
                <truncated>0</truncated>
                <difficult>0</difficult>
                <bndbox>
                        <xmin>{}</xmin>
                        <ymin>{}</ymin>
                        <xmax>{}</xmax>
                        <ymax>{}</ymax>
                </bndbox>
        </object>
'''

xml_footer = '''\
</annotation>
'''


def parse_annotaitons(fpath):
    
    img_dict = {}
    n = 0
    
    with open(fpath, 'r') as infh:
        for buf in infh:
            
            if n == 0:
                n = 1
                continue
            
            imageid, source, labelname, confidence, xmin, xmax, ymin, ymax, _, _, _, _, _ = buf.replace('\n', '').split(',')
             
            if imageid not in img_dict:
                img_dict[imageid] = []
            img_dict[imageid].append([float(xmin), float(xmax), float(ymin), float(ymax)])
    
    return img_dict



def download_image(fpath, img_dict):
   
    n = 0
    with open(fpath, 'r') as infh:
        for buf in infh:
            n = n + 1
            
            if n == 1:
                continue
            
            imagename, imageurl = buf.replace('\n', '').split(',')
            imageid = os.path.splitext(imagename)[0]
            if imageid in img_dict:
                print(imageurl)
                urllib.request.urlretrieve(imageurl, os.path.join('images', imagename))
                
                if n % 10 == 0:
                    time.sleep(5)
                

def make_labels(img_dict):
    
    for imageid, bboxes in img_dict.items():
        img_fname = os.path.join('images', imageid + '.jpg')
        txt_fname = os.path.join('images', imageid + '.xml')
        
        img = cv2.imread(img_fname)
        h, w, ch = img.shape
        
        with open(txt_fname, 'w') as outfh:
            _xml_header = xml_header.format(txt_fname, img_fname, w, h) + '\n'
            _xml_record = ''
            for xmin, xmax, ymin, ymax in bboxes:
                xmin = int(xmin * w)
                xmax = int(xmax * w)
                ymin = int(ymin * h)
                ymax = int(ymax * h)
                _xml_record = _xml_record + xml_record.format(xmin, ymin, xmax, ymax) + '\n'
           
            outfh.write(_xml_header + _xml_record + xml_footer)
    


if __name__ == '__main__':
    
    annotation_fpath = sys.argv[1]
    url_fpath = sys.argv[2]
    
    img_dict = parse_annotaitons(annotation_fpath)
    download_image(url_fpath, img_dict)
    make_labels(img_dict)          



