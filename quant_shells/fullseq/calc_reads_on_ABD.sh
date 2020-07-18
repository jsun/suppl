#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=2:00:00
#PBS -b 1
#PBS -m be
#PBS -N CALC_READS


PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment

cd ${PROJECT_DIR}
cd data


#
# HISAT stats
# 

log_fpath=n_reads_ABDgenome.tsv
rm ${log_fpath}

for bam_fpath in `ls bam | grep -v csi`
do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h bam/${bam_fpath} | grep 'chr.A' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    samtools view -h bam/${bam_fpath} | grep 'chr.B' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    samtools view -h bam/${bam_fpath} | grep 'chr.D' | grep 'NH:i:1' | wc -l >> ${log_fpath}
done


# 
# EAGLE-RC stats
# 

cd fullseq_eaglerc

log_fpath=${PROJECT_DIR}/data/fullseq_eaglerc/n_reads_ABDgenome.eaglerc.tsv
rm ${log_fpath}

cd eagle

for bam_dpath in `ls`
do
    echo ${bam_dpath} >> ${log_fpath}
    cd ${bam_dpath}
    
    for list_fpath in `ls *.ref.chr*.list`
    do
        echo -n "${bam_fpath}    " >> ${log_fpath}
         grep -c REF ${list_fpath} >> ${log_fpath}
    done
    
    cd ..
done




