# Dragonfly Classification ++



```bash
conda activate dragonflyenv

PROJECT_PATH=${HOME}/projects/dragonfly
PROJECTCLS_PATH=${PROJECT_PATH}/DragonflyPlus
DATA_PATH=${PROJECT_PATH}/data
SCRIPT_PATH=${PROJECT_PATH}/scripts
```



## Datasets

Check README.md in `data` to prepare datasets.




## Model trainin

```
Arch        dataset_W1   dataset_W2   datasset_F  dataset_W1F   dataset_W2F
mobilenet  
resnet     
vgg         
```


```bash
bash train_W1.sh
bash train_W2.sh
bash train_F.sh
bash train_W1F.sh
bash train_W2F.sh
```


























## Validation Summarise

To inference with test dataset, use the following samples.

```

cl=${DATA_PATH}/dataset_W1/class_labels.txt
val_dataset=/home/sonk414/projects/dragonfly/data/dataset_F2/images
val_results=./validation_results
model_names=("m1" "m2" "m3" "m4" "m5")

for model_name in ${model_names[@]}
do
    for dname in `ls ${val_dataset}`
    do
        # model i
        python predict.py -m ./weights/${model_name}_weights_run1.pth -c ${cl} -i ${val_dataset}/${dname} \
                          -o ./validation_results/${model_name}_weights_run1.${dname}.tsv
    done
done



cd validation_results
for model_name in ${model_names[@]}
do
    head -n 1 ${model_name}_weights_run1.Aeschnophlebia_longistigma.tsv > ${model_name}_weights_run1.all.tsv
    for fpath in `ls ${model_name}_weights_run1.* | grep -v all`
    do
        tail -n +2 ${fpath} >> ${model_name}_weights_run1.all.tsv
    done
done




# inference with single image
python predict.py -m ./weights/dragonfly-cls-model.v1.pth \
                  -i ./testimagedir/Sympetrum_darwinianum___4.cropped.png
```


