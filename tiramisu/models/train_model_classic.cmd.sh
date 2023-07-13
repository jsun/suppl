#!/bin/bash

# bash train_model_classic.cmd.sh 0 > qs_tiramisu_classic.log.0  2>&1


export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4
export VECLIB_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4



SGE_TASK_ID=$1

# project path
PROJECT_PATH=/data/ai_plantdisease/tiramisu
PROJECT_PATH=/home/sonk414/projects/tiramisu
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
            output_fpath=${result_dpath}/${alg}____${feature_type}____${random_type}____${disease_name}.tsv
            echo ${output_fpath}
            if [ -e ${output_fpath} ]; then
                echo "exists, skipped."
            else
                printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${feature_type}\e[m\n"
                ${PYTHON} model_classic.py    \
                    --algorithm ${alg}        \
                    --dataset ${disease_path} \
                    --feature-type ${feature_type}  \
                    --output ${output_fpath}
                printf "\e[31m# ^^^^^^^^^^ \e[m\n"
            fi
        done
    done
done



