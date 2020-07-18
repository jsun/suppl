#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=64gb
#PBS -l elapstim_req=72:00:00
#PBS -N WGENOME_MAP_HISAT

nCPU=16

PROJECT_DIR=/home/jqsun/research/wheat_genome
cd ${PROJECT_DIR}

BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc
DATA_DIR=${PROJECT_DIR}/data
FASTQ_DIR=${DATA_DIR}/fastq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bam
COUNTS_DIR=${DATA_DIR}/counts
GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2


mkdir -p ${BAM_DIR}
cd ${CLEAN_FASTQ_DIR}
for fastq_filepath in `ls *.clean.fastq.gz`; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`

    # mapping
    echo 'start mapping ...'
    echo `date`
    time ${BIN}/hisat2 -p ${nCPU} --no-spliced-alignment -x ${GENOME_INDEX} -U ${fastq_filepath} -S ${BAM_DIR}/${fastq_filename}.sam
    echo `date`
    echo 'finished.'

    # sorting
    cd ${BAM_DIR}
    ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
    ${BIN}/samtools index -c ${fastq_filename}.bam
        
    # remove sam
    rm ${fastq_filename}.sam
    cd -
done

