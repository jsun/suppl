#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -b 1
#PBS -m be
#PBS -N SHORTAGE_FULSEQ2


PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment
/home/jqsun/research/wheat_cold_treatment/data/clean_fullseq

cd ${PROJECT_DIR}
cd data/clean_fullseq


for fpath in `ls *_CS_*.fastq.gz | grep R1 | grep -v cold`
do
     ~/.pyenv/versions/ngs/bin/python ${PROJECT_DIR}/scripts/shortage_fastq.py ${fpath}
done


