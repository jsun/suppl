# Tiramisu


```bash
PROJECT_PATH=~/projects/tiramisu
```



# Data preparation


Put dataset into `formatted_data` and run the `format_data.R` scirpt.

```bash
cd ${PROJECT_PATH}
cd data

# create dataset for trainig classic models
for crop in cucumber eggplant strawberry tomato
do
    mkdir -p formatted_data/${crop}/shuffle
    python format_data.py "send210315/${crop}" "formatted_data/${crop}/shuffle" shuffle
    mkdir -p formatted_data/${crop}/randomize
    python format_data.py "send210315/${crop}" "formatted_data/${crop}/randomize" randomize
done

# create dataset for training deep nn mdoels


```




# Modeling

## Standard machine learning models

```bash
cd ${PROJECT_PATH}/models
mkdir test_results

for alg in svm rf dc knn lasso elasticnet
do
    for ft in category decimal
    do
        printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${ft}\e[m\n"

        python model_classic.py  \
        --algorithm ${alg} \
        --dataset ../data/send210315/cucumber/cucumber__udonkobyo__hompohasseimenseki.csv \
        --feature-type ${ft} \
        --output test_results/test_${alg}_${ft}.tsv \
        --test-run
    done
done

qsub train_model_classic.sh
```

## Deep neural network models

```bash
cd ${PROJECT_PATH}/models

dpath=${PROJECT_PATH}/data/send210315
mkdir -p ${dpath}/cv_dnnresults


# check cv training function for category type

python bake.py --mode cv --feature-type category --model L2 \
               --weight ${dpath}/cv_dnnresults/test_L2.cv_1.pth    \
               --train-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --valid-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --epochs 5 --batch-size 1024

python bake.py --mode cv --feature-type category --model L2 \
               --weight ${dpath}/cv_dnnresults/test_L2.cv_2.pth \
               --train-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --valid-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --epochs 5 --batch-size 1024


# check cv training function for decimal type

python bake.py --mode cv --feature-type decimal --model L2 \
               --weight ${dpath}/cv_dnnresults/test_L2.cv_1.pth    \
               --train-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --valid-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --epochs 5 --batch-size 1024

python bake.py --mode cv --feature-type decimal --model L2 \
               --weight ${dpath}/cv_dnnresults/test_L2.cv_2.pth \
               --train-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --valid-dataset ${dpath}/cucumber/cucumber__udonkobyo__hompohatsubyoyoritsu.csv \
               --epochs 5 --batch-size 1024

python summarise_cv.py L2 ${dpath}/cv_dnnresults ${dpath}/cv_dnnresults/test_bestparam.json.L2










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



# Summarization

```
R
source('scripts/eval_sunkishi_results.R')
```




