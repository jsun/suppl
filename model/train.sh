#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_cv
#$ -mods l_hard h_rt 144:00:00
#$ -t 1:3

source /etc/profile.d/modules.sh
module load cuda/10.0/10.0.130.1 cudnn/7.6.5


PROJECT_PATH=/home/sonk414/projects/tiramisu

cd ${PROJECT_PATH}

source ~/.profile
pyenv local tiramisu

cd model


model_arch=L${SGE_TASK_ID}

python train.py --model ${model_arch} --output ./weights/${model_arch}.pth       \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 1000 --batch-size 128 --train-mode cv



