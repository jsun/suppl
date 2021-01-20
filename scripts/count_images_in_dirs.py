import os
import sys
import glob

# python $this.py dragonfly_classes.txt /data/dataset_T/raw > dataset_T_summary.tsv

species_list_fpath = sys.argv[1]
indpath = sys.argv[2]

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
            n_images[os.path.basename(dpath)] += 1
    

for _species in species:
    print('{}\t{}'.format(_species.replace('_', ' '), n_images[_species]))



