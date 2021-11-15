#!/bin/bash

export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4
export VECLIB_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4


# bash train_model_dnn.cmd.sh 0 > qs_tiramisu_dnn.log.0  2>&1

SGE_TASK_ID=$1

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
            output_fpath=${result_dpath}/dnn____${ft}____${rt}____${disease_name}
            echo ${output_fpath}
            if [ -e ${output_fpath}_cv.10_relu_40-40_0.5.pth ]; then
                echo "exists, skipped"
            else
                ${PYTHON} ${PROJECT_PATH}/models/bake.py \
                    --mode cv --feature-type ${ft} --model L2 \
                    --weight ${result_dpath}/dnn____${ft}____${rt}____${disease_name}.pth \
                    --cv-dataset ${disease_path} \
                    --epochs 100 --batch-size 1024
            fi
        done
    done
done



