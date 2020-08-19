#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qs_tiramisu_cv
#$ -mods l_hard h_rt 144:00:00
#$ -t 1:30

# run mode
#   1: cross-validation to search best combination of parameters
#   2: train and test with the best combination of parameters
#   3: permutation test with the best combination of parameters
RUN_MODE=1

# crop name for analyses
CROP_NAME=Cucumber

# project path
PROJECT_PATH=~/projects/tiramisu




# load shell config data
# source /etc/profile.d/modules.sh
# module load cuda/10.0/10.0.130.1 cudnn/7.6.5
source ~/.profile




# move to project working space
cd ${PROJECT_PATH}
pyenv local tiramisu
cd model


# list up the target diseases
disease_names=($(ls -d ${PROJECT_PATH}/formatted_data/${CROP_NAME}__* | xargs -n1 basename))


# set the target disease in this thread
disease_name=${disease_names[${SGE_TASK_ID}]}
disease_path=${PROJECT_PATH}/formatted_data/${disease_name}
mkdir -p ${disease_path}/weights
mkdir -p ${disease_path}/weights_cv
mkdir -p ${disease_path}/weights_permutest_1
mkdir -p ${disease_path}/weights_permutest_2

best_param_fpath=${disease_path}/best_param.json

# cv to find the best parameters
if [ ${RUN_MODE} -eq 1 ]
then
    for model_arch in L1 L2
    do
        for dtrain_cv in `ls ${disease_path}/cv/train_*.tsv`
        do
            dvalid_cv=${dtrain_cv//train_/valid_}
            i=`basename ${dtrain_cv%.tsv}`
            i=${i//train_std_/}
            python bake.py --mode cv --model ${model_arch} \
                           --weight ${disease_path}/weights_cv/${model_arch}.cv_${i}.pth \
                           --train-dataset ${dtrain_cv} \
                           --valid-dataset ${dtrain_cv} 
        done
        
        # check the cv results and save best params into config file
        python summarise_cv.py ${model_arch} ${disease_path}/weights_cv ${best_param_fpath}.${model_arch}
    done
fi



if [ ${RUN_MODE} -eq 2 ]
then
    # some diseases may not contain training and validation data since there is no data before 2013,
    # if there is training in this disease, run for creating model
    if [ -e ${disease_path}/train_std.tsv ]
    then

    for model_arch in L1 L2
    do
        # train and validaiton with the best parameters
        python bake.py --mode train --model ${model_arch}                         \
                       --params ${best_param_fpath}.${model_arch}                 \
                       --weight ${disease_path}/weights/${model_arch}_best.pth    \
                       --train-dataset ${disease_path}/train_std.tsv              \
                       --valid-dataset ${disease_path}/valid_std.tsv
        
        python bake.py --mode valid --model ${model_arch}                            \
                       --params ${best_param_fpath}.${model_arch}                    \
                       --weight ${disease_path}/weights/${model_arch}_best.pth       \
                       --valid-dataset ${disease_path}/valid_std.tsv                 \
                       --output ${disease_path}/weights/${model_arch}_best.pth_valid
        
        # use all data to finalize model with the best parameters
        python bake.py ---mode train -model ${model_arch}                         \
                       --params ${best_param_fpath}.${model_arch}                 \
                       --weight ${disease_path}/weights/${model_arch}_final.pth   \
                       --train-dataset ${disease_path}/data_std.tsv               \
                       --valid-dataset ${disease_path}/data_std.tsv               \
        
        python bake.py --mode valid --model ${model_arch}                             \
                       --params ${best_param_fpath}.${model_arch}                     \
                       --weight ${disease_path}/weights/${model_arch}_final.pth       \
                       --valid-dataset ${disease_path}/data_std.tsv                   \
                       --output ${disease_path}/weights/${model_arch}_final.pth_valid
        
        python bake.py --mode valid --model ${model_arch}                            \
                       --params ${best_param_fpath}.${model_arch}                    \
                       --weight ${disease_path}/weights/${model_arch}_final.pth      \
                       --output ${disease_path}/weights/${model_arch}_final.pth_test \
                       --valid-dataset ${disease_path}/test_std.tsv
    done

    fi
fi



if [ ${RUN_MODE} -eq 3 ]
then
    for model_arch in L1 L2
    do
        for dtype in 1 2
        do
            for i in $(seq 1 100)
            do
                python bake.py --mode train --model ${model_arch} \
                               --params ${best_param_fpath}.${model_arch}                                  \
                               --weight ${disease_path}//weights_permutest_${dtype}/${model_arch}_${i}.pth \
                               --train-dataset ${disease_path}/randomize_${dtype}/train_${i}.tsv           \
                               --valid-dataset ${disease_path}/randomize_${dtype}/valid_${i}.tsv
        
                python bake.py --mode valid --model ${model_arch} \
                               --params ${best_param_fpath}.${model_arch}                                  \
                               --weight ${disease_path}/weights_permutest_${dtype}/${model_arch}_${i}.pth  \
                               --valid-dataset ${disease_path}//weights_permutest_${dtype}/valid_${i}.tsv  \
                               --output ${disease_path}//weights_permutest__${dtype}/${model_arch}_${i}.pth_valid
            done
        done
    done
fi





