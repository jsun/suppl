#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=64gb
#PBS -l elapstim_req=24:00:00
#PBS -N HB_TAG_DWSAMPL_QUANT
#PBS -t 2-8:2


nCPU=16
PROJECT_DIR=~/projects/HuoberBrezel
GENOME_DIR=~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.gff3


# 0.20, 0.40, 0.60, 0.80
# dw=0.${PBS_SUBREQNO}0
dwrate=('0.20' '0.40' '0.60' '0.80')


for dw in ${dwrate[@]}
do
    BATCH_DIR=${PROJECT_DIR}/data/tagseq_${dw}
    CLEAN_FASTQ_DIR=${BATCH_DIR}/clean_fastq
    BAM_DIR=${BATCH_DIR}/bam
    COUNTS_DIR=${BATCH_DIR}/counts

    cd ${PROJECT_DIR}
    mkdir -p ${BAM_DIR}
    
    cd ${CLEAN_FASTQ_DIR}
    for fastq_filepath in `ls *.clean.fastq.gz`; do
        cd ${CLEAN_FASTQ_DIR}
        echo ${fastq_filepath}
        fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`

        echo 'start mapping ...'
        echo `date`
        time hisat2 -p ${nCPU} --no-spliced-alignment -x ${GENOME_INDEX} \
                    -U ${fastq_filepath} -S ${BAM_DIR}/${fastq_filename}.sam
        echo `date`
        echo 'finished.'

        cd ${BAM_DIR}
        samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
        samtools index -c ${fastq_filename}.bam

        #rm ${fastq_filename}.sam
    done


    mkdir -p ${COUNTS_DIR}
    cd ${BAM_DIR}
    featureCounts -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                  -o ${COUNTS_DIR}/cs.counts.gene.ext_1.0k.tsv 20181109*.bam
    cd ${COUNTS_DIR}
    gzip *.tsv
    gzip *.tsv.summary
done




