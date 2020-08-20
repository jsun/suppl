#!/bin/bash

source ../local/global.sh


JOB_ROOT=${PROJECT_ROOT}/qc


FASTQ_REPORTS_BEFORE=${FASTQDIR}/reports
FASTQ_REPORTS_AFTER=${FASTQDIR}/reports_filtered


mkdir -p ${FASTQ_REPORTS_BEFORE}
mkdir -p ${FASTQ_REPORTS_AFTER}




cd ${FASTQDIR}

for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    ${SOURCEDIR}/FastQC/fastqc -t ${THREADS} --nogroup --noextract -o ${FASTQ_REPORTS_BEFORE} ${FASTQDIR}/${rnaseq_lib}_R1.fq.gz
    ${SOURCEDIR}/FastQC/fastqc -t ${THREADS} --nogroup --noextract -o ${FASTQ_REPORTS_BEFORE} ${FASTQDIR}/${rnaseq_lib}_R2.fq.gz
done



for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    java -Xmx4g -jar ${SOURCEDIR}/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE -threads ${THREADS} -phred33 \
        ${FASTQDIR}/${rnaseq_lib}_R1.fq.gz ${FASTQDIR}/${rnaseq_lib}_R2.fq.gz \
        ${FASTQDIR}/${rnaseq_lib}_filtered_R1.fq.gz \
        ${FASTQDIR}/${rnaseq_lib}_unpaired_R1.fq.gz \
        ${FASTQDIR}/${rnaseq_lib}_filtered_R2.fq.gz \
        ${FASTQDIR}/${rnaseq_lib}_unpaired_R2.fq.gz \
        ILLUMINACLIP:${SOURCEDIR}/Trimmomatic-0.36/adapters.fa:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done




for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    ${SOURCEDIR}/FastQC/fastqc -t ${THREADS} --nogroup --noextract -o ${FASTQ_REPORTS_AFTER} ${FASTQDIR}/${rnaseq_lib}_filtered_R1.fq.gz
    ${SOURCEDIR}/FastQC/fastqc -t ${THREADS} --nogroup --noextract -o ${FASTQ_REPORTS_AFTER} ${FASTQDIR}/${rnaseq_lib}_filtered_R2.fq.gz
done




rm ${FASTQ_REPORTS_BEFORE}/*.zip
rm ${FASTQ_REPORTS_AFTER}/*.zip


