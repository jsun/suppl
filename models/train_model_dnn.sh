#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_strawberry
#$ -mods l_hard h_rt 240:00:00
#$ -t 1:50


# crop name for analyses
# cucumber eggplant tomato strawberry
CROP_NAME=cucumber

# project path
PROJECT_PATH=~/projects/tiramisu
DATA_PATH=${PROJECT_PATH}/data/formatted_data
RESULT_PATH=${PROJECT_PATH}/data/cv_results

# load shell config data
#source ~/.profile

# move to project working space
cd ${PROJECT_PATH}
cd models
#pyenv local tiramisu



for rt in shuffle randomize
do
    data_path=${DATA_PATH}/${CROP_NAME}/${rt}4dnn
    disease_names=($(ls ${data_path}))
    disease_path=${data_path}/${disease_paths[${SGE_TASK_ID}]}
    disease_name=${disease_names[${SGE_TASK_ID}]}

    # path to save cv results
    result_dpath=${RESULT_PATH}/${CROP_NAME}/${rt}/dnn_models/${disease_name}
    mkdir -p ${result_dpath}
        
    # model architecture L1 and L2
    for ft in category decimal
    do
        for model_arch in L2
        do
            python ${PROJECT_PATH}/model/bake.py \
                --mode cv --feature-type ${rt} --model L2 \
                --weight ${result_dpath}/dnn____${ft}____${rt}____${disease_name}.pth \
                --cv-dataset ${disease_path} \
                --epochs 100 --batch-size 1024
        done
    done
done



