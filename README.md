# Tiramisu


```bash
PROJECT_PATH=~/projects/tiramisu
```



# Dataset

Put dataset into `formatted_data` and run the `format_data.classic.py` and `format_data.dnn.py`.

```bash
cd ${PROJECT_PATH}
cd data

# create dataset for trainig classic models

# extend ZIP from Kishi's dataset, and put all csv files into
# `formatted_data/random` and `formatted_data/shuffle` directories.
# format_data.classic.py is not used any more since Kishi shuffled and randomized dataset already.


# create dataset for training deep nn mdoels
cd formatted_data

mkdir -p random4dnn
for fpath in `ls random`
do
    mkdir -p "random4dnn/${fpath%_random.csv}"
    python ${PROJECT_PATH}/data/format_data.dnn.py "random/${fpath}" random4dnn/${fpath%_random.csv}
done
    
mkdir -p shuffle4dnn
for fpath in `ls shuffle`
do
    mkdir -p "shuffle4dnn/${fpath%_shuffle.csv}"
    python ${PROJECT_PATH}/data/format_data.dnn.py "shuffle/${fpath}" shuffle4dnn/${fpath%_shuffle.csv}
done

# python summary_datasets.py dataset_summary.xlsx
```




# Modeling

## Classic machine learning models

```bash
cd ${PROJECT_PATH}/models
mkdir -p test_results

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

# qsub is not working well for some deseases, no error-log, so cannot debug.
qsub train_model_classic.sh cucumber
qsub train_model_classic.sh eggplant
qsub train_model_classic.sh strawberry
qsub train_model_classic.sh tomato
```

## Deep neural network models

```bash
cd ${PROJECT_PATH}/models
dpath=${PROJECT_PATH}/data/formatted_data
mkdir -p test_results

# check cv training function for category type
python bake.py --mode cv --feature-type category --model L2 \
               --weight test_results/test_L2_category.pth   \
               --cv-dataset ${dpath}/cucumber/randomize4dnn/cucumber__udonkobyo__hompohatsubyoyoritsu \
               --epochs 5 --batch-size 1024

python bake.py --mode cv --feature-type decimal --model L2 \
               --weight test_results/test_L2_decimal.pth   \
               --cv-dataset ${dpath}/cucumber/randomize4dnn/cucumber__udonkobyo__hompohatsubyoyoritsu \
               --epochs 5 --batch-size 1024


qsub train_model_dnn.sh cucumber
qsub train_model_dnn.sh strawberry
qsub train_model_dnn.sh eggplant
qsub train_model_dnn.sh tomato
```


## Validation

Summarise validation results of DNN and classic models.

```
cd ${PROJECT_PATH}/models

python summarise_valid.py ${PROJECT_PATH}/data/cv_results/cucumber
python summarise_valid.py ${PROJECT_PATH}/data/cv_results/strawberry
python summarise_valid.py ${PROJECT_PATH}/data/cv_results/eggplant
python summarise_valid.py ${PROJECT_PATH}/data/cv_results/tomato

mkdir cv_results
cp ${PROJECT_PATH}/data/cv_results/cucumber/summary.tsv cv_results/cucumber_summary.tsv
cp ${PROJECT_PATH}/data/cv_results/strawberry/summary.tsv cv_results/strawberry_summary.tsv
cp ${PROJECT_PATH}/data/cv_results/eggplant/summary.tsv cv_results/eggplant_summary.tsv
cp ${PROJECT_PATH}/data/cv_results/tomato/summary.tsv cv_results/tomato_summary.tsv

R
> source('summarise_valid.R')
```


