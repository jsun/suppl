#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l elapstim_req=4:00:00
#PBS -b 1
#PBS -m be
#PBS -N NB_FUL_BAMSTATSTMP


PROJECT_DIR=~/projects/HuoberBrezel


cd ~/projects/HuoberBrezel/data
cd fullseq/bamstar

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
    grep 'chr.A' tmpbam_cdka8klxn > tmpbam_cdka8klxn.A
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.A >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.A >> ${log_fpath}
    
    echo '-- B --' >> ${log_fpath}
    grep 'chr.B' tmpbam_cdka8klxn > tmpbam_cdka8klxn.B
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.B >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.B >> ${log_fpath}
    
    echo '-- D --' >> ${log_fpath}
    grep 'chr.D' tmpbam_cdka8klxn > tmpbam_cdka8klxn.D
    samtools view -@ 8 -c -q 255 -f 3 tmpbam_cdka8klxn.D >> ${log_fpath}
    samtools view -@ 8 -c -q 255 -f 9 -F 4 tmpbam_cdka8klxn.D >> ${log_fpath}
    
    rm tmpbam_cdka8klxn tmpbam_cdka8klxn.A tmpbam_cdka8klxn.B tmpbam_cdka8klxn.D
done


