#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_dwsample
#$ -mods l_hard h_rt 240:00:00
#$ -t 2-8:2


PROJECT_DIR=~/projects/HuoberBrezel
CLEAN_FASTQ_PATH=${PROJECT_DIR}/data/tagseq/clean_fastq
SCRIPTS=${PROJECT_DIR}/quant_shells/scripts

cd ${PROJECT_DIR}

# 0.20, 0.40, 0.60, 0.80
dw=0.${SGE_TASK_ID}0


#for dw in ${dwrate[@]}
#do
    batch_dpath=${PROJECT_DIR}/data/tagseq_${dw}/clean_fastq
    mkdir -p ${batch_dpath}
    
    for fpath in `ls ${CLEAN_FASTQ_PATH}/*.gz | grep RS1`
    do
        fname=`basename ${fpath}`
        python ${SCRIPTS}/downsampling_fastq.py ${fpath} ${batch_dpath}/${fname} ${dw}
    done
    cd

#done

