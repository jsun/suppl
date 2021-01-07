import os
import sys
import re
import glob
import shutil
from PIL import Image
from PIL import ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True



def parse_xml(file_path):
        objects = []
        img_shape = [None, None, None]
        object_dict = {'name': None, 'xmin': None, 'ymin': None, 'xmax': None, 'ymax': None}
        
        # regex for searching object's information
        re_width = re.compile(r'<width>([0-9]+)</width>')
        re_height = re.compile(r'<height>([0-9]+)</height>')
        re_depth = re.compile(r'<depth>([0-9]+)</depth>')
        re_name = re.compile(r'<name>(.+)</name>')
        re_xmin = re.compile(r'<xmin>([0-9]+)</xmin>')
        re_ymin = re.compile(r'<ymin>([0-9]+)</ymin>')
        re_xmax = re.compile(r'<xmax>([0-9]+)</xmax>')
        re_ymax = re.compile(r'<ymax>([0-9]+)</ymax>')

        # open Pascal VOC XML and read object information
        with open(file_path, 'r') as xmlfh:

            is_object_record = False

            for line in xmlfh:
                if not is_object_record:
                    m = re_width.search(line)
                    if m:
                        img_shape[0] = int(m.group(1))
                    m = re_height.search(line)
                    if m:
                        img_shape[1] = int(m.group(1))
                    m = re_depth.search(line)
                    if m:
                        img_shape[2] = int(m.group(1))

                if '<object>' in line:
                    is_object_record = True
                    continue

                if '</object>' in line:
                    objects.append(object_dict)
                    object_dict = {'name': None, 'xmin': None, 'ymin': None, 'xmax': None, 'ymax': None}
                    is_object_record = False
                    continue

                if is_object_record:
                    m = re_name.search(line)
                    if m:
                        object_dict['name'] = m.group(1)

                    m = re_xmin.search(line)
                    if m:
                        object_dict['xmin'] = int(m.group(1))

                    m = re_ymin.search(line)
                    if m:
                        object_dict['ymin'] = int(m.group(1))

                    m = re_xmax.search(line)
                    if m:
                        object_dict['xmax'] = int(m.group(1))

                    m = re_ymax.search(line)
                    if m:
                        object_dict['ymax'] = int(m.group(1))

        return {'shape': tuple(img_shape), 'objects': objects}



def crop_images(img_fpath):
    
    objects = parse_xml(os.path.splitext(img_fpath)[0] + '.xml')
    im = Image.open(img_fpath)
    for obj in objects['objects']:
        obj = (obj['xmin'], obj['ymin'], obj['xmax'], obj['ymax'])
        im_cropped = im.crop(obj)
        output_fpath = os.path.splitext(img_fpath)[0] + '__cropped__' + '_'.join(map(str, obj)) + '.jpg'
        
        if im.info.get('exif') is not None:
            im_cropped.save(output_fpath, quality=100, exif=im.info.get('exif'))





def crop(input_dirpath):
    
    for image_path in sorted(glob.glob(os.path.join(input_dirpath, '*', '*.jpg'))):
        if not os.path.exists(os.path.splitext(image_path)[0] + '.xml'):
            continue
        
        crop_images(image_path)
        


   


if __name__ == '__main__':
    
    input_dpath = sys.argv[1]
    crop(input_dpath)
        


