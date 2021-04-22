#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_classic_cucumber
#$ -mods l_hard h_rt 240:00:00
#$ -t 1:50

# run mode
RUN_MODE=1

# crop name for analyses
# cucumber eggplant tomato strawberry
CROP_NAME=cucumber

# project path
PROJECT_PATH=~/projects/tiramisu
DATA_PATH=${PROJECT_PATH}/data/send210315

# load shell config data
#source ~/.profile

# move to project working space
cd ${PROJECT_PATH}
pyenv local tiramisu


# list up the target diseases
disease_names=($(ls -d ${DATA_PATH}/${CROP_NAME}/${CROP_NAME}__* | xargs -n1 basename))

# set the target disease in this thread
disease_datapath=${DATA_PATH}/${CROP_NAME}/${disease_names[${SGE_TASK_ID}]}
disease_name=${disease_names[${SGE_TASK_ID}]%.csv}

# cv results
result_dpath=${DATA_PATH}/cv_results/${CROP_NAME}
mkdir -p ${result_dpath}


# cv to find the best parameters
if [ ${RUN_MODE} -eq 1 ]
then
    
    for alg in svm rf dc knn lasso elasticnet
    do
        for ft in category decimal
        do
            printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${ft}\e[m\n"

            python model_classic.py  \
                --algorithm ${alg} \
                --dataset ${disease_datapath} \
                --feature-type ${ft} \
                --output ${result_dpath}/${alg}____${ft}____${disease_name}.tsv
        done
    done

fi




