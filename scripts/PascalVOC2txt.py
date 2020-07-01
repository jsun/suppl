import sys
import os
import xml.etree.ElementTree as ET



def pascalVOC2txt(xml_fpath):

    fdata = []

    with open(xml_fpath, 'r') as xmlfh:
        
        xml_root = ET.parse(xmlfh).getroot()
        img_path = xml_root.find('path').text
        
        for obj in xml_root.iter('object'):
            
            xml_bbox = obj.find('bndbox')
            _fdata = [img_path,
                     xml_bbox.find('xmin').text, xml_bbox.find('ymin').text, xml_bbox.find('xmax').text, xml_bbox.find('ymax').text,
                     'dragonfly']
            fdata.append(','.join(_fdata))
    return fdata




if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        raise ValueError('pytho $this image_dpath output_fpath')
    
    image_dpath = sys.argv[1]
    output_fpath = sys.argv[2]
    
    fdata = []
    
    for r, d, f in os.walk(image_dpath):
        for _f in f:
            if _f.endswith('.xml'):
                fdata.extend(pascalVOC2txt(os.path.join(r, _f)))
    
    with open(output_fpath, 'w') as outfh:
        outfh.write('\n'.join(fdata))



