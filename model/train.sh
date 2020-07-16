#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_cv
#$ -mods l_hard h_rt 144:00:00
#$ -t 1:3

# 1: cross-validation search
# 2: train and test with default nn
# 3: permutation test with default nn
task_id=1

source /etc/profile.d/modules.sh
module load cuda/10.0/10.0.130.1 cudnn/7.6.5


PROJECT_PATH=~/projects/tiramisu

cd ${PROJECT_PATH}

source ~/.profile
pyenv local tiramisu

cd model


model_arch=L${SGE_TASK_ID}

if [ ${task_id} -eq 1 ]
then
    python bake.py --model ${model_arch} --weight ./weights/${model_arch}.pth       \
               --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
               --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
               --epochs 1000 --batch-size 128 --mode cv
fi


if [ ${task_id} -eq 2 ]
then
    python bake.py --model L2 --weight ./weights/jppnet.pth \
               --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
               --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
               --epochs 1000 --batch-size 128 --mode train

    python bake.py --model L2 --weight ./weights/jppnet.pth --output valid_outputs/validresult     \
               --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
               --mode valid
    
    python bake.py --model L2 --weight ./weights/jppnet.pth --output valid_outputs/testresult     \
               --valid-dataset ../formatted_data/kyuuri_honpo_percent.test.tsv \
               --mode valid
    
fi



if [ ${task_id} -eq 3 ]
then
    
    for dtype in Type1 Type2
    do
    
        mkdir -p permutest/${dtype}
        
        for i in $(seq 1 1000)
        do
            fpath=${PROJECT_PATH}/formatted_data/randomized_${dtype}/kyuuri_honpo_percent.train.${dtype}.${i}.tsv

            python bake.py --weight ./permutest/${dtype}.pth \
                   --train-dataset ${fpath} \
                   --valid-dataset ${PROJECT_PATH}/formatted_data/kyuuri_honpo_percent.valid.tsv \
                   --epochs 1000 --batch-size 128 --mode train
        
            python bake.py --weight ./permutest/${dtype}.pth --output permutest/${dtype}/p${i}       \
                   --valid-dataset ${PROJECT_PATH}/formatted_data/kyuuri_honpo_percent.valid.tsv \
                   --mode valid
        
        done


    done



fi





