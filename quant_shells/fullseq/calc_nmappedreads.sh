#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=4:00:00
#PBS -b 1
#PBS -m be
#PBS -N NB_FUL_BAMSTATS


PROJECT_DIR=~/projects/HuoberBrezel


cd ~/projects/HuoberBrezel/data
cd fullseq/bam

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'chr.A' | grep -c 'NH:i:1' >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'chr.B' | grep -c 'NH:i:1' >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'chr.D' | grep -c 'NH:i:1' >> ${log_fpath}
done



cd ~/projects/HuoberBrezel/data
cd fullseq/bamstar

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -q 255 -f 3 ${bam_fpath} | grep 'chr.A' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 -f 3 ${bam_fpath} | grep 'chr.B' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 -f 3 ${bam_fpath} | grep 'chr.D' | grep -v "@SQ" | wc -l >> ${log_fpath}
done






cd ~/projects/HuoberBrezel/data
cd fullseq/eaglerc

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrA.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrB.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrD.list >> ${log_fpath}
done



