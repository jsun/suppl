# Dragonfly Classification


```bash
PROJECT_PATH=${HOME}/projects/dragonfly
PROJECTCLS_PATH=${PROJECT_PATH}/DragonflyCls
DATA_PATH=${PROJECT_PATH}/data
SCRIPT_PATH=${PROJECT_PATH}/scripts

export OMP_NUM_THREADS=1 # avoid warnings in macOS
```



## Datasets

Check README.md in `data` to prepare datasets.



## Model training

The 6 archtectures x 5 types of dataset.

```bash
python train.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                --model-arch vgg \
                --model-outpath ./weights_species/F__vgg__test.pth \
                --traindata ./test_datasets/W1 \
                --validdata ./test_datasets/W1 \
                -e 3 -b 32 -l 0.001

```

```bash
bash train_W1.sh
bash train_W2.sh
bash train_F.sh
bash train_W1F.sh
bash train_W2F.sh
```

Then find the best model and rename to `best.pth`.


```bash
ln weights_species/xxx.pth weights_species/best_resnet152.pth
```


## Inference

```bash
# image foler
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/best.pth \
                  -i testimagedir \
                  -o test_result.tmp

# single image
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/best.pth \
                  -i testimagedir/Sympetrum_darwinianum_1.png \
                  -o test_result.tmp

```


# DragonflyCls + DragonflyMesh


```bash
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/best.pth \
                  --mesh ${DATA_PATH}/meshmatrix.tsv.gz \
                  -i testimagedir/Sympetrum_darwinianum_2.jpg \
                  -o test_result.tmp


python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/best.pth \
                  --mesh ${DATA_PATH}/meshmatrix.tsv.gz \
                  -i testimagedir \
                  -o test_result.tmp
```



# Validation with field photo

for species level

```bash
cd ${PROJECTCLS_PATH}
dpath=${DATA_PATH}/dataset_T/cropped_image

for model in `ls weights_species/*.pth`
do
    for d in `ls ${dpath}`
    do
        echo "${model%.pth} -- ${dpath}/${d}"
        model_arch=(${model//__/ })
        model_arch=${model_arch[1]}
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.image.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 50         \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d50.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 100           \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d100.tsv  \
                      --overwrite
    done
done

```


for genus level

```bash
cd ${PROJECTCLS_PATH}
dpath=${DATA_PATH}/dataset_Tg/cropped_image

for model in `ls weights_genus/*.pth`
do
    for d in `ls ${dpath}`
    do
        echo "${model%.pth} -- ${dpath}/${d}"
        model_arch=(${model//__/ })
        model_arch=${model_arch[1]}
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.image.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 50           \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d50.tsv     \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 100          \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d100.tsv  \
                      --overwrite
    done
done
```




