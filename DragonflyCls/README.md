# Dragonfly Classification


```bash
PROJECT_PATH=${HOME}/projects/dragonfly
PROJECTCLS_PATH=${PROJECT_PATH}/DragonflyCls
DATA_PATH=${PROJECT_PATH}/data
SCRIPT_PATH=${PROJECT_PATH}/scripts
```



## Datasets

Check README.md in `data` to prepare datasets.



## Model training

The 6 archtectures x 5 types of dataset.


```bash
bash train_W1.sh
bash train_W2.sh
bash train_F.sh
bash train_W1F.sh
bash train_W2F.sh
```

Then find the best model and rename to `best.pth`.


```bash
mv weights/xxx.pth weights/best.pth
```


## Inference

```bash
python predict.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                  --model-arch resnet152 \
                  --model-weight weights/best.pth \
                  -i testimagedir/Sympetrum_darwinianum___4.jpg \
                  -o test_result.tmp

```


