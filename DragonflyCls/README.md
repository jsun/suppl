# DragonflyCls

To create deep learning models for odonates prediction, run the following scripts one by one.



## Setup Variables

```bash
PROJECT_PATH=${HOME}/projects/dragonfly
PROJECTCLS_PATH=${PROJECT_PATH}/DragonflyCls
DATA_PATH=${PROJECT_PATH}/data
SCRIPT_PATH=${PROJECT_PATH}/scripts
```



## Dataset Preparation

Check README.md in `data` to prepare datasets.



## Training

To train VGG16 architecture with dataset W1, run the following scripts.

```bash
python train.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                --model-arch vgg \
                --model-outpath ./weights_species/F__vgg__test.pth \
                --traindata ./test_datasets/W1 \
                --validdata ./test_datasets/W1 \
                -e 3 -b 32 -l 0.001

```

Run shell scripts to train six architectures (VGG16, VGG19, ResNet18, ResNet152, DenseNet, MobileNet)
with the five types of dataset (W1, W2, W1F, W2F, F).


```bash
bash train_W1.sh
bash train_W2.sh
bash train_F.sh
bash train_W1F.sh
bash train_W2F.sh
```


## Inference

To perform odonate prediciton with trained models, use the following scripts.

```bash
cd ${PROJECTCLS_PATH}

# image foler
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/W1__resnet152__1.pth \
                  -i testimagedir \
                  -o test_result.tmp

# single image
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/W1__resnet152__1.pth \
                  -i testimagedir/Sympetrum_darwinianum_1.png \
                  -o test_result.tmp
```

To perform odonate prediction with trained image models and occurrence records,
use the following scripts.

```bash
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/W1__resnet152__1.pth \
                  --mesh ${DATA_PATH}/meshmatrix.tsv.gz \
                  -i testimagedir/Sympetrum_darwinianum_2.jpg \
                  -o test_result.tmp


python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights_species/W1__resnet152__1.pth \
                  --mesh ${DATA_PATH}/meshmatrix.tsv.gz \
                  -i testimagedir \
                  -o test_result.tmp
```



## Validation

To validate trained models with/without occurrence records with dataset T,
run the following scripts.

### Species Level 

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
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 50          \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d50.tsv    \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 100         \
                      -i ${dpath}/${d} -o ${model%.pth}.T_valid.d100.tsv   \
                      --overwrite
    done
done

```

### Genus Level

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
                      -i ${dpath}/${d} -o ${model%.pth}.Tg_valid.image.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 50           \
                      -i ${dpath}/${d} -o ${model%.pth}.Tg_valid.d50.tsv     \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 100          \
                      -i ${dpath}/${d} -o ${model%.pth}.Tg_valid.d100.tsv  \
                      --overwrite
    done
done
```


### Summarization

Use R script to summarize training status and validation results.

```bash
cd ${PROJECTCLS_PATH}
Rscript eval_results.R
```




## Validation (for Specimen dataset)

To validate trained models with/without occurrence records with dataset T,
run the following scripts.

### Species Level 

```bash
cd ${PROJECTCLS_PATH}
dpath=${DATA_PATH}/dataset_Tw/raw

for model in `ls weights_species/*.pth`
do
    for d in `ls ${dpath}`
    do
        echo "${model%.pth} -- ${dpath}/${d}"
        model_arch=(${model//__/ })
        model_arch=${model_arch[1]}
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      -i ${dpath}/${d} -o ${model%.pth}.Tw_valid.image.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 50          \
                      -i ${dpath}/${d} -o ${model%.pth}.Tw_valid.d50.tsv    \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      --mesh ${DATA_PATH}/meshmatrix.tsv.gz -d 100         \
                      -i ${dpath}/${d} -o ${model%.pth}.Tw_valid.d100.tsv   \
                      --overwrite
    done
done

```

### Genus Level

```bash
cd ${PROJECTCLS_PATH}
dpath=${DATA_PATH}/dataset_Twg/raw

for model in `ls weights_genus/*.pth`
do
    for d in `ls ${dpath}`
    do
        echo "${model%.pth} -- ${dpath}/${d}"
        model_arch=(${model//__/ })
        model_arch=${model_arch[1]}
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}   \
                      -i ${dpath}/${d} -o ${model%.pth}.Twg_valid.image.tsv  \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 50           \
                      -i ${dpath}/${d} -o ${model%.pth}.Twg_valid.d50.tsv     \
                      --overwrite
        
        python predict.py --class-label ${DATA_PATH}/dragonflyg_classes.txt \
                      --model-arch ${model_arch} --model-weight ${model}    \
                      --mesh ${DATA_PATH}/meshmatrixg.tsv.gz -d 100          \
                      -i ${dpath}/${d} -o ${model%.pth}.Twg_valid.d100.tsv  \
                      --overwrite
    done
done
```




