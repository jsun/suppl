#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_classic
#$ -mods l_hard h_rt 720:00:00
#$ -t 1-210:1


# project path
PROJECT_PATH=/data/ai_plantdisease/tiramisu
DATA_PATH=${PROJECT_PATH}/data/formatted_data
RESULT_PATH=${PROJECT_PATH}/data/cv_results
PYTHON=/home/sonk414/.pyenv/versions/tiramisu/bin/python

# move to project working space
cd ${PROJECT_PATH}
cd models


# randomized or shuffled data
for random_type in shuffle random
do
    data_dpath=${DATA_PATH}/${random_type}
    
    # list up the target diseases
    disease_paths=($(ls -d ${data_dpath}/*.csv))
    disease_names=($(ls -d ${data_dpath}/*.csv | xargs -n1 basename))
    
    # set the target disease in this thread
    disease_path=${disease_paths[${SGE_TASK_ID}]}
    disease_name=${disease_names[${SGE_TASK_ID}]}
    disease_name=${disease_name%.csv}
    disease_name=${disease_name%_random}
    disease_name=${disease_name%_shuffle}

    # path to save cv results
    result_dpath=${RESULT_PATH}/classic_models
    mkdir -p ${result_dpath}
    
    # cv to find the best parameters
    echo "${disease_name}"
    
    for alg in svm rf dc knn lasso elasticnet
    do
        for feature_type in category decimal
        do
                printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${feature_type}\e[m\n"
                ${PYTHON} model_classic.py    \
                    --algorithm ${alg}        \
                    --dataset ${disease_path} \
                    --feature-type ${feature_type}  \
                    --output ${result_dpath}/${alg}____${feature_type}____${random_type}____${disease_name}.tsv
                printf "\e[31m# ^^^^^^^^^^ \e[m\n"
            done
        done
done



