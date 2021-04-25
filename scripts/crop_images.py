import os
import sys
import re
import glob
import shutil
from PIL import Image
from PIL import ImageFile
from PIL import ExifTags
import PIL.ImageOps


ImageFile.LOAD_TRUNCATED_IMAGES = True


def get_jpeg_info(img_fpath):
    lat = None
    lng = None
    capture_date = None
    im = Image.open(img_fpath)

    exif = im._getexif()

    if exif is not None:
        exif = {ExifTags.TAGS[k]: v for k, v in exif.items() if k in ExifTags.TAGS}
        if 'GPSInfo' in exif:
            gps_tags = exif['GPSInfo']
            gps = {ExifTags.GPSTAGS.get(t, t): gps_tags[t] for t in gps_tags}
            is_lat = 'GPSLatitude' in gps
            is_lat_ref = 'GPSLatitudeRef' in gps
            is_lon = 'GPSLongitude' in gps
            is_lon_ref = 'GPSLongitudeRef' in gps

            if is_lat and is_lat_ref and is_lon and is_lon_ref:
                lat = gps['GPSLatitude']
                lat_ref = gps['GPSLatitudeRef']
                if lat_ref == 'N':
                    lat_sign = 1.0
                elif lat_ref == 'S':
                    lat_sign = -1.0
                lon = gps['GPSLongitude']
                lon_ref = gps['GPSLongitudeRef']
                if lon_ref == 'E':
                    lon_sign = 1.0
                elif lon_ref == 'W':
                    lon_sign = -1.0
                lat = lat_sign * lat[0] + lat[1] / 60 + lat[2] / 3600
                lng = lon_sign * lon[0] + lon[1] / 60 + lon[2] / 3600

        if 'DateTimeOriginal' in exif:
            capture_date = exif['DateTimeOriginal']
            capture_date = capture_date.split(' ')[0].replace(':', '-')
    return (capture_date, lat, lng)




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



def crop_images(img_fpath, mode):
    objects = parse_xml(os.path.splitext(img_fpath)[0] + '.xml')
    im = Image.open(img_fpath)
    im = PIL.ImageOps.exif_transpose(im)
    for obj in objects['objects']:
        obj = (obj['xmin'], obj['ymin'], obj['xmax'], obj['ymax'])
        im_cropped = im.crop(obj)
        output_fpath = os.path.splitext(img_fpath)[0] + '__cropped__' + '_'.join(map(str, obj)) + '.jpg'
        if mode == 'exif':
            im_cropped.save(output_fpath, quality=100, exif=im.info.get('exif'))
        else:
            im_cropped.save(output_fpath, quality=100)


def crop(mode, input_dirpath):
    for image_path in sorted(glob.glob(os.path.join(input_dirpath, '*', '*.jpg'))):
        if not os.path.exists(os.path.splitext(image_path)[0] + '.xml'):
            continue
        
        if mode == 'exif':
            exif_info = get_jpeg_info(image_path)
            print(image_path, exif_info)
            if (exif_info[0] is None) or (exif_info[1] is None) or (exif_info[2] is None):
                print('No EXIF: ' + image_path)
                continue
        
        crop_images(image_path, mode)
        


   


if __name__ == '__main__':
    
    mode = sys.argv[1]
    input_dpath = sys.argv[2]
    crop(mode, input_dpath)
        


