import os
import sys
import glob
import shutil
from pandabox.imgUtils import imgUtils
from PIL import ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True



def crop(input_dirpath, output_dirpath):
    
    iu = imgUtils()
    
    for image_path in sorted(glob.glob(os.path.join(input_dirpath, '*', '*.jpg'))):
        
        _output_dirpath = os.path.join(output_dirpath, os.path.basename(os.path.dirname(image_path)))
        if not os.path.exists(_output_dirpath):
            os.makedirs(_output_dirpath)
        
        try:
            xml_path = image_path.replace('.jpg', '.xml')
            obj = iu.parse_PascalVOC(xml_path)
            iu.crop_images(image_path, obj, _output_dirpath)
        
            for dname in os.listdir(_output_dirpath):
                if os.path.isdir(os.path.join(_output_dirpath, dname)):
                    for f in glob.glob(os.path.join(_output_dirpath, dname, '*.png')):
                        print(image_path)
                        shutil.move(f, os.path.join(os.path.dirname(os.path.dirname(f)), os.path.basename(f)))
                    os.rmdir(os.path.dirname(f))
        except:
            pass



   


if __name__ == '__main__':
    
    input_dpath = sys.argv[1]
    output_dpath = sys.argv[2]
    crop(input_dpath, output_dpath)
        


      


