#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_train_W1F
#$ -mods l_hard h_rt 144:00:00
#$ -t 1:10

source /etc/profile.d/modules.sh
module load cuda/10.0/10.0.130.1 cudnn/7.6.5


PROJECT_PATH=${HOME}/projects/dragonfly
PROJECTCLS_PATH=${PROJECT_PATH}/DragonflyCls
DATA_PATH=${PROJECT_PATH}/data
SCRIPT_PATH=${PROJECT_PATH}/scripts



cd ${PROJECTCLS_PATH}


source ~/.profile
pyenv local dragonfly


i=${SGE_TASK_ID}


model_archs=(vgg resnet mobilenet vgg19 resnet152 densenet)

for model_arch in "${model_archs[@]}"
do
    echo ${model_arch}

    time python train.py --class-label ${DATA_PATH}/dragonfly_classes.txt \
                --model-arch ${model_arch} \
                --model-outpath ./weights/W1F__${model_arch}__${i}.pth \
                --traindata ${DATA_PATH}/dataset_W1F/merged \
                --validdata ${DATA_PATH}/dataset_T/cropped_image  \
                -e 50 -b 32 -l 0.001
done

