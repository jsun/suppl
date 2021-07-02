#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_dnn
#$ -mods l_hard h_rt 720:00:00
#$ -t 1-90:1


# crop name for analyses
# cucumber eggplant tomato strawberry
CROP_NAME=$1

# project path
# PROJECT_PATH=~/projects/tiramisu
PROJECT_PATH=/data/ai_plantdisease/tiramisu
DATA_PATH=${PROJECT_PATH}/data/formatted_data
RESULT_PATH=${PROJECT_PATH}/data/cv_results
PYTHON=/data/ai_plantdisease/tiramisu/python_env/bin/python


# move to project working space
cd ${PROJECT_PATH}
cd models



for rt in shuffle randomize
do
    i=`expr ${SGE_TASK_ID} - 1`

    data_path=${DATA_PATH}/${CROP_NAME}/${rt}4dnn
    disease_names=($(ls ${data_path}))
    disease_path=${data_path}/${disease_names[${i}]}
    disease_name=${disease_names[${i}]}

    # path to save cv results
    result_dpath=${RESULT_PATH}/${CROP_NAME}/${rt}/dnn_models/${disease_name}
    mkdir -p ${result_dpath}
        
    # model architecture L1 and L2
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



