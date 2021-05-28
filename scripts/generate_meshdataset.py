import os
import sys
import re
import numpy as np
import pandas as pd


level = sys.argv[1]
kankyo_fpath = sys.argv[2]
spname_fpath = sys.argv[3]
class_fpath = sys.argv[4]
output_fpath = sys.argv[5]


def mesh2gps(mesh_code):
    mesh_code = str(mesh_code)
    lat = int(mesh_code[0:2]) * 2 / 3
    lng = int(mesh_code[2:4]) + 100
    
    if len(mesh_code) > 4:
        if len(mesh_code) >= 6:
            lat += int(mesh_code[4]) * 2 / 3 / 8
            lng += int(mesh_code[5]) / 8

    return (lat, lng)




# get class labels (the order should be matched to image-model outputs)
class_labels = []
with open(class_fpath, 'r') as infh:
    for buf in infh:
        class_labels.append(buf.replace('\n', ''))


# get metadata to convert species ID to species biname
id2class = {}
with open(spname_fpath, 'r') as infh:
    infh.readline()
    for buf in infh:
        bufs = buf.replace('\n', '').split(',')
        id2class[bufs[0]] = bufs[4] + '_' + bufs[6]

# read Kankyosho public data
# and manually modifiy Kankyosho data according to rearrangement of taxonomic orders
## Rhipidolestes okinawanus: 392722, 392746, 392756, 392757, 392860, 392870
## Rhipidolestes shozoi:     392860, 392870, 402801, 402811, 402812
## Rhipidolestes amamiensis: 412857, 412867, 422922, 222932, 422933, 422944, 473002
## Rhipidolestes asatoi:     472935, 472945
## Anotogaster klossi:       362336, 362337, 362346, 362347, 362441, 362451
## Rhipidolestes yakusimensis: remove 472935, 472945, and add 473002 from the original set
## Anotogaster sieboldii:      remove 362336, 362337, 362346, 362347, 362441, 362451 from the original set
fdata_mesh = []
fdata_species = []
with open(kankyo_fpath, 'r') as infh:
    for buf in infh:
        bufs = buf.replace('\n', '').split(',')
        cl = id2class[bufs[0]]
        if cl == 'Rhipidolestes_yakusimensis':
            if bufs[1] in  ['472935', '472945']:
                print('removed: ' + cl + ' -- ' + bufs[1])
        elif cl == 'Anotogaster_sieboldii':
            if bufs[1] in ['362336', '362337', '362346', '362347', '362441', '362451']:
                print('removed: ' + cl + ' -- ' + bufs[1])
        else:
            fdata_mesh.append(bufs[1])
            fdata_species.append(id2class[bufs[0]])

fdata_species.extend(['Rhipidolestes_okinawanus'] * 6)
fdata_mesh.extend(['392722', '392746', '392756', '392757', '392860', '392870'])
fdata_species.extend(['Rhipidolestes_shozoi'] * 5)
fdata_mesh.extend(['392860', '392870', '402801', '402811', '402812'])
fdata_species.extend(['Rhipidolestes_amamiensis'] * 7)
fdata_mesh.extend(['412857', '412867', '422922', '222932', '422933', '422944', '473002'])
fdata_species.extend(['Rhipidolestes_asatoi'] * 2)
fdata_mesh.extend(['472935', '472945'])
fdata_species.extend(['Anotogaster_klossi'] * 6)
fdata_mesh.extend(['362336', '362337', '362346', '362347', '362441', '362451'])
fdata_species.extend(['Rhipidolestes_yakusimensis'])
fdata_mesh.extend(['473002'])



# change species name (level) to genus name (level)
if level == 'genus':
    for i, spname in enumerate(fdata_species):
        fdata_species[i] = spname.split('_')[0]


# mesh to lat&lng
latlng = []
for _fdata_mesh in sorted(list(set(fdata_mesh))):
    latlng.append(mesh2gps(_fdata_mesh))

latlng = pd.DataFrame(latlng, columns=['lat', 'lng'],
                      index=sorted(list(set(fdata_mesh))))


# make appearance matrix
print(len(class_labels))
dmat = pd.DataFrame(np.zeros((len(set(fdata_mesh)), len(class_labels))))
dmat.columns = class_labels
dmat.index = sorted(list(set(fdata_mesh)))
# appearance matrix summary
dsum = pd.DataFrame(np.zeros((len(set(fdata_mesh)), len(class_labels))))
dsum.columns = class_labels
dsum.index = sorted(list(set(fdata_mesh)))

for _mesh, _species in zip(fdata_mesh, fdata_species):
    if _species in class_labels:
        dmat.loc[_mesh, _species] = 1
        dsum.loc[_mesh, _species] += 1

dmat = pd.concat([latlng, dmat], axis=1)


dsum = dsum.sum(axis=0)
print(dsum)

# write out the data
dmat.to_csv(output_fpath, header=True, index=True, sep='\t', compression='gzip')
dsum.to_csv(output_fpath.replace('.tsv', '').replace('.gz', '') + '.summary.tsv', header=False, index=True, sep='\t')



