# Data preparation

```
# R
source('scripts/format_data.R')
```




# Modeling


```
cd nn_model

# test
python bake.py --model L1 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --mode train

python bake.py --model L1 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --mode cv

python bake.py --model L2 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --mode cv

python bake.py --model L3 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --mode cv


# grid search
qsub train.sh

```








