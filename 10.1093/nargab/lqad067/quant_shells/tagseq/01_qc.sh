#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -m be
#PBS -N HB_TAG_QC


nCPU=8
PROJECT_DIR=~/projects/HuoberBrezel
DATA_DIR=${PROJECT_DIR}/data
FASTQ_DIR=${PROJECT_DIR}/data/tagseq/fastq
CLEAN_FASTQ_DIR=${PROJECT_DIR}/data/tagseq/clean_fastq


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
    echo ${fastq_path}
    
    java -jar trimmomatic-0.39.jar  \
              SE -phred33 -threads ${nCPU} \
              ${fastq_path} \
              ${CLEAN_FASTQ_DIR}/${fastq_name}.clean.trimmomatic.fastq.gz          \
              ILLUMINACLIP:/tools/local/meta/trimmomatic_adapters.txt:2:30:10 \
              TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:20 MINLEN:40
done


# poly-A trimming
cd ${CLEAN_FASTQ_DIR}
python ../../../quant_shells/scripts/trim_polyA.py . . clean.trimmomatic.fastq.gz clean.fastq.gz ${nCPU}
# rm *.clean.trimmomatic.fastq.gz


# check the qualities of the original files
cd ${CLEAN_FASTQ_DIR}
mkdir -p fastqc_reports

fastqc --nogroup -t ${nCPU}                     \
       -a ~/local/meta/trimmomatic_adapters.txt \
       -o fastqc_reports                        \
       --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -



