#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -b 1
#PBS -m be
#PBS -N WCT_DWSAMPL_FASTQ



PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment/data
cd ${PROJECT_DIR}


dwdata=(0.20 0.40 0.60 0.80)
dwdata=(0.20)

for dw in ${dwdata[@]}
do
    batch_dpath=clean_tagseq_${dw}
    mkdir -p ${batch_dpath}/fastq
    cd ${batch_dpath}/fastq
    for fpath in `ls ../../clean_tagseq_cs/*.gz`
    do
        output_fpath=`basename ${fpath}`
        echo $fpath
        echo $output_fpath
        python ../../../scripts/downsampling_fastq.py ${fpath} ${output_fpath} ${dw}
    done
    cd

done

