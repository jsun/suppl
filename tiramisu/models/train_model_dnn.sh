#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_dnn
#$ -mods l_hard h_rt 720:00:00
#$ -t 1-210:1


# project path
PROJECT_PATH=/data/ai_plantdisease/tiramisu
PROJECT_PATH=/home/sonk414/projects/tiramisu
DATA_PATH=${PROJECT_PATH}/data/formatted_data
RESULT_PATH=${PROJECT_PATH}/data/cv_results
PYTHON=/home/sonk414/.pyenv/versions/tiramisu/bin/python


cd ${PROJECT_PATH}
cd models


for rt in shuffle random
do
    data_path=${DATA_PATH}/${rt}4dnn
   
    disease_paths=($(ls -d ${data_path}/*))
    disease_names=($(ls -d ${data_path}/* | xargs -n1 basename))

    disease_path=${disease_paths[${SGE_TASK_ID}]}
    disease_name=${disease_names[${SGE_TASK_ID}]}
    
    # path to save cv results
    result_dpath=${RESULT_PATH}/dnn_models/${rt}/${disease_name}
    mkdir -p ${result_dpath}
        
    for ft in category decimal
    do
        for model_arch in L2
        do
            ${PYTHON} ${PROJECT_PATH}/models/bake.py \
                --mode cv --feature-type ${ft} --model L2 \
                --weight ${result_dpath}/dnn____${ft}____${rt}____${disease_name}.pth \
                --cv-dataset ${disease_path} \
                --epochs 100 --batch-size 1024
        done
    done
done



