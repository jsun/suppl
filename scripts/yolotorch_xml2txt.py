import sys
import os
import glob
from pandabox import imgUtils

def _xml2txt(input_fpath, output_fpath):
    imgutils = imgUtils()
    xmlinfo = imgutils.parse_PascalVOC(input_fpath)
    img_w, img_h, img_ch = xmlinfo['shape']
    with open(output_fpath, 'w') as outfh:
        for obj in xmlinfo['objects']:
            class_id = None
            obj_x_center = None
            obj_y_center = None
            obj_width = None
            obj_height = None

            if obj['name'] == 'dragonfly':
                class_id = 0
            else:
                raise ValueError('found unexpect label: ' + input_fpath)

            xmin = obj['xmin'] / img_w
            xmax = obj['xmax'] / img_w
            ymin = obj['ymin'] / img_h
            ymax = obj['ymax'] / img_h
            obj_x_center = (xmin + xmax) / 2
            obj_y_center = (ymin + ymax) / 2
            obj_width = xmax - xmin
            obj_height = ymax - ymin

            dat_record = '{} {} {} {} {}\n'.format(class_id,
                        obj_x_center, obj_y_center, obj_width, obj_height)
            outfh.write(dat_record)


def xml2txt(input_dpath):
    for xml_path in glob.glob(input_dpath + '/*.xml'):
        txt_path = os.path.splitext(xml_path)[0] + '.txt'
        _xml2txt(xml_path, txt_path)



if __name__ == '__main__':
    input_dpath = sys.argv[1]
    xml2txt(input_dpath)


