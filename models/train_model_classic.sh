#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_classic
#$ -mods l_hard h_rt 720:00:00
#$ -t 1-90:1


# set crop name: cucumber eggplant tomato strawberry
# CROP_NAME=strawberry
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




# randomized or shuffled data
for rt in shuffle randomize
do
    data_path=${DATA_PATH}/${CROP_NAME}/${rt}

    # list up the target diseases
    disease_paths=($(ls -d ${data_path}/*.csv))
    disease_names=($(ls -d ${data_path}/*.csv | xargs -n1 basename))

    # set the target disease in this thread
    disease_path=${disease_paths[${SGE_TASK_ID}]}
    disease_name=${disease_names[${SGE_TASK_ID}]}
    disease_name=${disease_name%.csv}
    disease_name=${disease_name%.randomize}
    disease_name=${disease_name%.shuffle}

    # path to save cv results
    result_dpath=${RESULT_PATH}/${CROP_NAME}/${rt}/classic_models
    mkdir -p ${result_dpath}

    # cv to find the best parameters
    echo "${disease_name}"
    echo "${disease_datapath}"
    
    for alg in svm rf dc knn lasso elasticnet
    do
        echo "#####"
        for ft in category decimal
        do
                printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${ft}\e[m\n"
                ${PYTHON} model_classic.py  \
                    --algorithm ${alg}        \
                    --dataset ${disease_path} \
                    --feature-type ${ft}      \
                    --randomize-type ${rt}    \
                    --output ${result_dpath}/${alg}____${ft}____${rt}____${disease_name}.tsv
            done
        done
done



