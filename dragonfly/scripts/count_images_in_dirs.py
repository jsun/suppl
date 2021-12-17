import os
import sys
import glob

from PIL import Image
from PIL import ExifTags
from PIL.ExifTags import TAGS


# python $this.py dragonfly_classes.txt /data/dataset_T/raw > dataset_T_summary.tsv


species_list_fpath = sys.argv[1]
indpath = sys.argv[2]






def get_gps(fpath):
    im = Image.open(fpath)
    lat = None
    lng = None
    try:
        exif = im._getexif()
        exif = {
            ExifTags.TAGS[k]: v for k, v in exif.items() if k in ExifTags.TAGS
        }
        if 'GPSInfo' in exif:
            gps_tags = exif['GPSInfo']
            gps = {
                ExifTags.GPSTAGS.get(t, t): gps_tags[t] for t in gps_tags
            }
            is_lat = 'GPSLatitude' in gps
            is_lat_ref = 'GPSLatitudeRef' in gps
            is_lng = 'GPSLongitude' in gps
            is_lng_ref = 'GPSLongitudeRef' in gps
            if is_lat and is_lat_ref and is_lng and is_lng_ref:
                lat = gps['GPSLatitude']
                lat_ref = gps['GPSLatitudeRef']
                if lat_ref == 'N':
                    lat_sign = 1.0
                elif lat_ref == 'S':
                    lat_sign = -1.0
                lng = gps['GPSLongitude']
                lng_ref = gps['GPSLongitudeRef']
                if lng_ref == 'E':
                    lng_sign = 1.0
                elif lng_ref == 'W':
                    lng_sign = -1.0
                lat = lat_sign * lat[0] + lat[1] / 60 + lat[2] / 3600
                lng = lng_sign * lng[0] + lng[1] / 60 + lng[2] / 3600
    except:
        pass
        # print(fpath + ' has no EXIF!')

    return lat, lng



def get_captured_datetime(fpath):
    captured_datetime = None
    im = Image.open(fpath)

    exif = im._getexif()
    try:
        for exif_id, exif_val in exif.items():
            tag = TAGS.get(exif_id, exif_id)
            if tag == 'DateTimeOriginal':
                captured_datetime = exif_val
    except:
        pass
        # print(fpath + ' has no EXIF!')

    im.close()

    return captured_datetime







species = []
with open(species_list_fpath, 'r') as infh:
    for line in infh:
        species.append(line.replace('\n', ''))


n_images = {}
for _species in species:
    n_images[_species] = 0


for dpath in sorted(glob.glob(os.path.join(indpath, '*'))):
    n = 0
    for fpath in glob.glob(os.path.join(dpath, '*')):
        if os.path.splitext(fpath)[1].lower() in ['.jpg', '.jpeg', '.png']:
            if 'dataset_T' in dpath:
                if get_captured_datetime(fpath) is not None and get_gps(fpath) is not None:
                    n_images[os.path.basename(dpath)] += 1
            else:
                n_images[os.path.basename(dpath)] += 1
    

for _species in species:
    print('{}\t{}'.format(_species.replace('_', ' '), n_images[_species]))


print('#' + str(sum(n_images.values())))

