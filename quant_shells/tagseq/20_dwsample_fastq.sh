#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -b 1
#PBS -m be
#PBS -N HB_TAG_DWSAMPL
#PBS -t 2-8:2



PROJECT_DIR=~/projects/HuoberBrezel
CLEAN_FASTQ_PATH=${PROJECT_DIR}/data/tagseq/clean_fastq
SCRIPTS=${PROJECT_DIR}/quant_shells/scripts

cd ${PROJECT_DIR}


# 0.20, 0.40, 0.60, 0.80
dwrate=0.${PBS_SUBREQNO}0


for dw in ${dwrate[@]}
do
    batch_dpath=${PROJECT_DIR}/data/tagseq_${dw}/clean_fastq
    mkdir -p ${batch_dpath}
    
    for fpath in `ls ${CLEAN_FASTQ_PATH}/*.gz | grep RS1`
    do
        fname=`basename ${fpath}`
        python ${SCRIPTS}/downsampling_fastq.py ${fpath} ${batch_dpath}/${fname} ${dw}
    done
    cd

done

