#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_full_bamstats
#$ -mods l_hard h_rt 288:00:00



PROJECT_DIR=~/projects/HuoberBrezel


cd ${PROJECT_DIR}/data/fullseq/bam

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'Chr.A' | grep -c 'NH:i:1' >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'Chr.B' | grep -c 'NH:i:1' >> ${log_fpath}
    samtools view -h -f 0x2 ${bam_fpath} | grep 'Chr.D' | grep -c 'NH:i:1' >> ${log_fpath}
done







cd ${PROJECT_DIR}/data/shortfullseq/bam

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h ${bam_fpath} | grep 'Chr.A' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    samtools view -h ${bam_fpath} | grep 'Chr.B' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    samtools view -h ${bam_fpath} | grep 'Chr.D' | grep 'NH:i:1' | wc -l >> ${log_fpath}
done






cd ${PROJECT_DIR}/data/fullseq/bamstar

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -@ 8 ${bam_fpath} > tmpbam_cdka8klxn

    echo '-- A + B + D + U --' >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn >> ${log_fpath}

    echo '-- A --' >> ${log_fpath}
    grep 'Chr.A' tmpbam_cdka8klxn > tmpbam_cdka8klxn.A
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.A >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.A >> ${log_fpath}

    echo '-- B --' >> ${log_fpath}
    grep 'Chr.B' tmpbam_cdka8klxn > tmpbam_cdka8klxn.B
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.B >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.B >> ${log_fpath}

    echo '-- D --' >> ${log_fpath}
    grep 'Chr.D' tmpbam_cdka8klxn > tmpbam_cdka8klxn.D
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.D >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.D >> ${log_fpath}

    rm tmpbam_cdka8klxn tmpbam_cdka8klxn.A tmpbam_cdka8klxn.B tmpbam_cdka8klxn.D
done






cd ${PROJECT_DIR}/data/fullseq/eaglerc

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for list_dpath in `ls`; do
    echo ${list_dpath} >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrA.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrB.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrD.list >> ${log_fpath}
done



