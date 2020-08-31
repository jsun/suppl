#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -m be
#PBS -N HB_FUL_QC


pQC1=0
pTrimLQ=1
pQC2=0



PROJECT_DIR=~/projects/HuoberBrezel
cd ${PROJECT_DIR}


nCPU=8
BIN=~/local/bin
UTILS=~/local/utilfunc


DATA_DIR=${PROJECT_DIR}/data
FASTQ_DIR=${PROJECT_DIR}/data/fullseq/fastq
CLEAN_FASTQ_DIR=${PROJECT_DIR}/data/fullseq/clean_fastq






# check the qualities of the original files
if [ ${pQC1} -eq 1 ]; then
    cd ${FASTQ_DIR}
    mkdir -p fastqc_reports
    # python ${UTILS}/parse_fastq.py .
    ${BIN}/fastqc --nogroup -t ${nCPU}                     \
                  -a ~/local/meta/trimmomatic_adapters.txt \
                  -o fastqc_reports                        \
                  --kmers 7 *_R1.fastq.gz
    rm fastqc_reports/*.zip
    cd -
fi


# quality filtering and poly-A trimming
if [ ${pTrimLQ} -eq 1 ]; then
    cd ${DATA_DIR}
    mkdir -p ${CLEAN_FASTQ_DIR}

    for fastq_path in `ls ${FASTQ_DIR}/*_R1.fastq.gz | grep W2017_1_CS-4`; do
        fastq_name=`basename ${fastq_path} _R1.fastq.gz`
        fastq1_path=${FASTQ_DIR}/${fastq_name}_R1.fastq.gz
        fastq2_path=${FASTQ_DIR}/${fastq_name}_R2.fastq.gz
        
        java -jar ${BIN}/trimmomatic-0.36.jar   \
              PE -phred33 -threads ${nCPU}  \
              ${fastq1_path} ${fastq2_path} \
              ${FASTQ_DIR}/${fastq_name}_R1.cleaned.fastq.gz  \
              ${FASTQ_DIR}/${fastq_name}_R1.unpaired.fastq.gz \
              ${FASTQ_DIR}/${fastq_name}_R2.cleaned.fastq.gz  \
              ${FASTQ_DIR}/${fastq_name}_R2.unpaired.fastq.gz \
              ILLUMINACLIP:/home/jqsun/local/meta/trimmomatic_adapters.txt:2:30:10 \
              TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:20 MINLEN:40
    done
    
    mv ${FASTQ_DIR}/*_cleaned.fastq.gz ${CLEAN_FASTQ_DIR}/
fi




# check the qualities of the original files
if [ ${pQC2} -eq 1 ]; then
    cd ${CLEAN_FASTQ_DIR}
    mkdir -p fastqc_reports
    # python ${UTILS}/parse_fastq.py .
    ${BIN}/fastqc --nogroup -t ${nCPU}                     \
                  -a ~/local/meta/trimmomatic_adapters.txt \
                  -o fastqc_reports                        \
                  --kmers 7 *.fastq.gz
    rm fastqc_reports/*.zip
    cd -
fi




