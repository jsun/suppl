#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -m be
#PBS -N HB_FUL_QC



nCPU=8
PROJECT_DIR=~/projects/HuoberBrezel
DATA_DIR=${PROJECT_DIR}/data
FASTQ_DIR=${PROJECT_DIR}/data/fullseq/fastq
CLEAN_FASTQ_DIR=${PROJECT_DIR}/data/fullseq/clean_fastq


cd ${PROJECT_DIR}


# check the qualities of the original files
cd ${FASTQ_DIR}
mkdir -p fastqc_reports
fastqc --nogroup -t ${nCPU}                     \
       -a ~/local/meta/trimmomatic_adapters.txt \
       -o fastqc_reports                        \
       --kmers 7 *_R1.fastq.gz
rm fastqc_reports/*.zip
cd -


# quality filtering and poly-A trimming
cd ${DATA_DIR}
mkdir -p ${CLEAN_FASTQ_DIR}

for fastq_path in `ls ${FASTQ_DIR}/*_R1.fastq.gz`; do
    fastq_name=`basename ${fastq_path} _R1.fastq.gz`
    fastq1_path=${FASTQ_DIR}/${fastq_name}_R1.fastq.gz
    fastq2_path=${FASTQ_DIR}/${fastq_name}_R2.fastq.gz
     
    java -jar trimmomatic-0.39.jar   \
              PE -phred33 -threads ${nCPU}  \
              ${fastq1_path} ${fastq2_path} \
              ${FASTQ_DIR}/${fastq_name}_R1.cleaned.fastq.gz  \
              ${FASTQ_DIR}/${fastq_name}_R1.unpaired.fastq.gz \
              ${FASTQ_DIR}/${fastq_name}_R2.cleaned.fastq.gz  \
              ${FASTQ_DIR}/${fastq_name}_R2.unpaired.fastq.gz \
              ILLUMINACLIP:/tools/local/meta/trimmomatic_adapters.txt:2:30:10 \
              TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:20 MINLEN:40
done
mv ${FASTQ_DIR}/*.cleaned.fastq.gz ${CLEAN_FASTQ_DIR}/


# check the qualities of the original files
cd ${CLEAN_FASTQ_DIR}
mkdir -p fastqc_reports
fastqc --nogroup -t ${nCPU}                     \
       -a ~/local/meta/trimmomatic_adapters.txt \
       -o fastqc_reports                        \
       --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -

