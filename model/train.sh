#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_eggplant
#$ -mods l_hard h_rt 240:00:00
#$ -t 1:50

# run mode
RUN_MODE=1

# crop name for analyses
# cucumber eggplant tomato strawberry
CROP_NAME=eggplant

# project path
PROJECT_PATH=~/projects/tiramisu




# load shell config data
source ~/.profile




# move to project working space
cd ${PROJECT_PATH}
pyenv local tiramisu


# list up the target diseases
disease_names=($(ls -d ${PROJECT_PATH}/formatted_data/${CROP_NAME}__* | xargs -n1 basename))


# set the target disease in this thread
disease_name=${disease_names[${SGE_TASK_ID}]}
disease_path=${PROJECT_PATH}/formatted_data/${disease_name}



# cv to find the best parameters
if [ ${RUN_MODE} -eq 1 ]
then
    cd ${disease_path}
    
    # model architecture L1 and L2
    for model_arch in L2
    do
        for dataset in norm_datasets null_datasets
        do
            cd ${disease_path}/${dataset}
            
            # from 001 to 100
            for i in `ls`
            do
                cd ${disease_path}/${dataset}/${i}
                
                # cross-validation for fixing hyper-parameters
                # --------------------------------------------
                cd cv
                for dtrain_cv in `ls train_*.tsv`
                do
                    dvalid_cv=${dtrain_cv//train_/valid_}
                    di=`basename ${dtrain_cv%.tsv}`
                    di=${di//train_std_/}
                    
                    python ${PROJECT_PATH}/model/bake.py --mode cv --model ${model_arch} \
                            --weight weight_${model_arch}.cv_${di}.pth \
                            --train-dataset ${dtrain_cv} \
                            --valid-dataset ${dvalid_cv} 
                done
                cd ..
                # check the cv results and save best params into config file
                python ${PROJECT_PATH}/model/summarise_cv.py ${model_arch} cv param.${model_arch}
                
                # train and validation
                # --------------------
                python ${PROJECT_PATH}/model/bake.py --mode train --model ${model_arch} \
                            --params param.${model_arch}          \
                            --weight weight_${model_arch}.pth     \
                            --train-dataset train_std.tsv         \
                            --valid-dataset valid_std.tsv
                
                python ${PROJECT_PATH}/model/bake.py --mode valid --model ${model_arch} \
                            --params param.${model_arch}          \
                            --weight weight_${model_arch}.pth     \
                            --valid-dataset valid_std.tsv         \
                            --output modelvalid_${model_arch}
                
            done  # dataset i (001, 002, ..., 100) 
        done # dataset (norm, null)
    done # model (L1, L2)
fi




# summarise the results
mkdir -p results/${CROP_NAME}
mkdir -p results/${CROP_NAME}/fig

rm results/${CROP_NAME}/valid_summary.tsv
cd ${PROJECT_PATH}/model
for disease_name in ${disease_names[@]}
do
    echo ${disease_name}
    python summarise_valid.py ${PROJECT_PATH}/formatted_data/${disease_name} >> results/${CROP_NAME}/valid_summary.tsv
    cp ${PROJECT_PATH}/formatted_data/${disease_name}/rmse_hist.png results/${CROP_NAME}/fig/${disease_name}.hist.png
done




