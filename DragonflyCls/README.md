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
                --model-outpath ./weights/F__vgg__test.pth \
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
mv weights/xxx.pth weights/best_resnet152.pth
```


## Inference

```bash
# image foler
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights/best.pth \
                  -i testimagedir \
                  -o test_result.tmp

# single image
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights/best.pth \
                  -i testimagedir/Sympetrum_darwinianum_1.png \
                  -o test_result.tmp
```


# DragonflyCls + DragonflyMesh


```bash
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights/best.pth \
                  --mesh ${DATA_PATH}/meshmatrix.tsv.gz \
                  -i testimagedir/Sympetrum_darwinianum___4.jpg \
                  -o test_result.tmp

```







