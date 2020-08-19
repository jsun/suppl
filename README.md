# Tiramisu


```bash
PROJECT_PATH=~/projects/tiramisu
```



# Data preparation


Put dataset into `formatted_data` and run the `format_data.R` scirpt.

```bash
cd ${PROJECT_PATH}
Rscript scripts/format_data.R 
```




# Modeling

Check the scripts run well.


```bash
cd ${PROJECT_PATH}
cd model

dpath=${PROJECT_PATH}/formatted_data/test
mkdir -p ${dpath}/weights
mkdir -p ${dpath}/weights_cv


# check cv training function
python bake.py --mode cv --model L1 \
               --weight ${dpath}/weights_cv/L1.cv_1.pth \
               --train-dataset ${dpath}/cv/train_std_1.tsv \
               --valid-dataset ${dpath}/cv/valid_std_1.tsv \
               --epochs 5 --batch-size 1024

python bake.py --mode cv --model L1 \
               --weight ${dpath}/weights_cv/L1.cv_2.pth \
               --train-dataset ${dpath}/cv/train_std_2.tsv \
               --valid-dataset ${dpath}/cv/valid_std_2.tsv \
               --epochs 5 --batch-size 1024

python summarise_cv.py L1 ${dpath}/weights_cv ${dpath}/best_param.json.L1




# check trainning function
python bake.py --mode train --model L1 \
               --params ${dpath}/best_param.json.L1 \
               --weight ${dpath}/weights/L1_best.pth \
               --train-dataset ${dpath}/train_std.tsv    \
               --valid-dataset ${dpath}/valid_std.tsv    \
               --epochs 5 --batch-size 1024

python bake.py --mode valid --model L1 \
               --params ${dpath}/best_param.json.L1 \
               --weight ${dpath}/weights/L1_best.pth \
               --valid-dataset ${dpath}/valid_std.tsv    \
               --output ${dpath}/weights/L1_best.pth_valid \


# grid search
qsub train.sh

```








