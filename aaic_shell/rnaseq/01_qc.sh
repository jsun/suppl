#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -m be
#PBS -N WCT_QC




BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc

PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment
DATA_DIR=${PROJECT_DIR}/data
FULLSEQ_DIR=${DATA_DIR}/fullseq
TAQSEQ_DIR=${DATA_DIR}/tagseq
CLEAN_FULLSEQ_DIR=${DATA_DIR}/clean_fullseq
CLEAN_TAQSEQ_DIR=${DATA_DIR}/clean_tagseq
nCPU=8

cd ${PROJECT_DIR}




# check the qualities of the original files
cd ${FULLSEQ_DIR}
mkdir fastqc_reports
python ${UTILS}/parse_fastq.py .
${BIN}/fastqc --nogroup -t ${nCPU} \
              -a ~/local/meta/trimmomatic_adapters.txt \
              -o fastqc_reports \
              --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -
cd ${TAQSEQ_DIR}
mkdir fastqc_reports
python ${UTILS}/parse_fastq.py .
${BIN}/fastqc --nogroup -t ${nCPU} \
              -a ~/local/meta/trimmomatic_adapters.txt \
              -o fastqc_reports \
              --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -




# clearning reads (fullseq)
cd ${PROJECT_DIR}
cd ${DATA_DIR}
mkdir -p ${CLEAN_FULLSEQ_DIR}

for fastq_path in `ls ${FULLSEQ_DIR}/*_R1.fastq.gz`; do
    fastq_name=`basename ${fastq_path} _R1.fastq.gz`
    fastq1_path=${FULLSEQ_DIR}/${fastq_name}_R1.fastq.gz
    fastq2_path=${FULLSEQ_DIR}/${fastq_name}_R2.fastq.gz
    echo ${fastq_path}

    java -jar ${BIN}/trimmomatic-0.36.jar   \
              PE -phred33 -threads ${nCPU}  \
              ${fastq1_path} ${fastq2_path} \
              ${CLEAN_FULLSEQ_DIR}/${fastq_name}_R1.cleaned.fastq.gz  \
              ${CLEAN_FULLSEQ_DIR}/${fastq_name}_R1.unpaired.fastq.gz \
              ${CLEAN_FULLSEQ_DIR}/${fastq_name}_R2.cleaned.fastq.gz  \
              ${CLEAN_FULLSEQ_DIR}/${fastq_name}_R2.unpaired.fastq.gz \
              ILLUMINACLIP:/home/jqsun/local/meta/trimmomatic_adapters.txt:2:30:10 \
              TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:20 MINLEN:40
done




cd ${PROJECT_DIR}
cd ${DATA_DIR}
mkdir -p ${CLEAN_TAQSEQ_DIR}

for fastq_path in `ls ${TAQSEQ_DIR}/*_R1.fastq.gz`; do
    fastq_name=`basename ${fastq_path} _R1.fastq.gz`
    echo ${fastq_path}

    java -jar ${BIN}/trimmomatic-0.36.jar  \
              SE -phred33 -threads ${nCPU} \
              ${fastq_path} \
              ${CLEAN_TAQSEQ_DIR}/${fastq_name}.clean.trimmomatic.fastq.gz         \
              ILLUMINACLIP:/home/jqsun/local/meta/trimmomatic_adapters.txt:2:30:10 \
              TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:20 MINLEN:40

    python ${UTILS}/trim_polyA.py ${CLEAN_TAQSEQ_DIR}/${fastq_name}.clean.trimmomatic.fastq.gz \
                                  ${CLEAN_TAQSEQ_DIR}/${fastq_name}.cleaned.fastq.gz
    rm ${CLEAN_TAQSEQ_DIR}/${fastq_name}.clean.trimmomatic.fastq.gz
done









# check the qualities of the original files
cd ${CLEAN_FULLSEQ_DIR}
mkdir fastqc_reports
${BIN}/fastqc --nogroup -t ${nCPU} \
              -a ~/local/meta/trimmomatic_adapters.txt \
              -o fastqc_reports_beforeQC \
              --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -
cd ${CLEAN_TAQSEQ_DIR}
mkdir fastqc_reports
${BIN}/fastqc --nogroup -t ${nCPU} \
              -a ~/local/meta/trimmomatic_adapters.txt \
              -o fastqc_reports_beforeQC \
              --kmers 7 *.fastq.gz
rm fastqc_reports/*.zip
cd -




