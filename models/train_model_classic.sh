#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_classic_cucumber
#$ -mods l_hard h_rt 240:00:00
#$ -t 1:50


# set crop name: cucumber eggplant tomato strawberry
CROP_NAME=cucumber

# project path
PROJECT_PATH=~/projects/tiramisu
DATA_PATH=${PROJECT_PATH}/data/formatted_data
RESULT_PATH=${PROJECT_PATH}/data/cv_results

# load shell config data
#source ~/.profile

# move to project working space
cd ${PROJECT_PATH}
#pyenv local tiramisu
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
    disease_name=${disease_names[${SGE_TASK_ID}]%.csv}
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
                ~/.pyenv/versions/tiramisu/bin/python model_classic.py  \
                    --algorithm ${alg}        \
                    --dataset ${disease_path} \
                    --feature-type ${ft}      \
                    --randomize-type ${rt}    \
                    --output ${result_dpath}/${alg}____${ft}____${rt}____${disease_name}.tsv
            done
        done
done



