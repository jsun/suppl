#!/bin/bash


PROJECT_DIR=~/projects/HuoberBrezel


 
cd ~/projects/HuoberBrezel/data
for dpath in `ls -d tagseq*`; do
    cd ${dpath}/bam

    log_fpath=nreads_mapped_stats.tsv
    if [ -f ${log_fpath} ]; then
        rm ${log_fpath}
    fi

    for bam_fpath in `ls *.bam`; do
        echo ${bam_fpath} >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.A' | grep 'NH:i:1' | wc -l >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.B' | grep 'NH:i:1' | wc -l >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.D' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    done

    cd -
done




cd ~/projects/HuoberBrezel/data
cd tagseq/bamstar

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.A' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.B' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.D' | grep -v "@SQ" | wc -l >> ${log_fpath}
done





cd ~/projects/HuoberBrezel/data
cd tagseq/eaglerc

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for list_dpath in `ls -d *`; do
    echo ${list_dpath} >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrA.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrB.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrD.list | grep -c REF >> ${log_fpath}
done


cd tagseq/eaglercngi

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for list_dpath in `ls -d *`; do
    echo ${list_dpath} >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrA.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrB.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrD.list | grep -c REF >> ${log_fpath}
done






