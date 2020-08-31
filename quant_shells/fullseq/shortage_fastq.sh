#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -b 1
#PBS -m be
#PBS -N HB_FUL_SHORTFQ


PROJECT_DIR=~/projects/HuoberBrezel


cd ${PROJECT_DIR}
cd data/fullseq/clean_fastq


for fpath in `ls *.fastq.gz | grep R1`
do
     ~/.pyenv/versions/ngs/bin/python ${PROJECT_DIR}/quant_shells/scripts/shortage_fastq.py ${fpath}
done




