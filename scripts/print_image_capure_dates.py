import os
import sys
import glob
from PIL import Image
from PIL import ExifTags
from PIL.ExifTags import TAGS
 

# pwd
# # ~/projects/dragonfly/data/dataset_T
# python ../../scripts/print_image_capure_dates.py cropped_image
# python ../../scripts/print_image_capure_dates.py raw
 

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
 

def main(target_path):
    n_images = 0
    n_noexif = 0
    datetimes = []
    for fpath in glob.glob(os.path.join(target_path, '*', '*.jpg')):
        n_images += 1
        datetime = get_captured_datetime(fpath)
        gps = get_gps(fpath)
        if gps[0] is not None:
            if datetime is not None:
                datetimes.append(datetime)
            print('{}\t{}\t{}\t{}'.format(fpath, datetime, gps[0], gps[1]))
        else:
            n_noexif += 1
        
        
    print('#{}-{}'.format(sorted(datetimes)[0], sorted(datetimes)[-1]))
    print('#{}, {}, {}, {}'.format(n_images, n_noexif, n_images - n_noexif))



 
if __name__ == '__main__':
    
    target_path = sys.argv[1]
    main(target_path)


