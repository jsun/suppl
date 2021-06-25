#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_full_shortagereads
#$ -mods l_hard h_rt 12:00:00



PROJECT_DIR=~/projects/HuoberBrezel


cd ${PROJECT_DIR}
cd data/fullseq/clean_fastq


for fpath in `ls *.fastq.gz | grep R1`
do
     ~/.pyenv/versions/bio/bin/python ${PROJECT_DIR}/quant_shells/scripts/shortage_fastq.py ${fpath}
done


