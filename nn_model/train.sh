#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_jppnet_cv
#$ -mods l_hard h_rt 72:00:00

source /etc/profile.d/modules.sh
module load cuda/10.0/10.0.130.1 cudnn/7.6.5


PROJECT_PATH=/home/sonk414/projects/jppnet

cd ${PROJECT_PATH}

source ~/.profile
pyenv local jppnet

cd nn_model


python train.py --model L1 --output ./weights/L1.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 200 --batch-size 1024 --train-mode cv

python train.py --model L2 --output ./weights/L2.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 200 --batch-size 1024 --train-mode cv

python train.py --model L3 --output ./weights/L3.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 200 --batch-size 1024 --train-mode cv



