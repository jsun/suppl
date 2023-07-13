#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_dwsamplequant
#$ -mods l_hard h_rt 240:00:00
#$ -t 2-8:2



nCPU=16
PROJECT_DIR=~/projects/HuoberBrezel
BIN=~/local/bin


# 0.20, 0.40, 0.60, 0.80
dw=0.${SGE_TASK_ID}0
BATCH_DIR=${PROJECT_DIR}/data/tagseq_${dw}

CLEAN_FASTQ_DIR=${BATCH_DIR}/clean_fastq
BAM_DIR=${BATCH_DIR}/bam

GENOME_DIR=~/projects/db/genome/IWGSC_RefSeq_v2.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
GTF_FILEPATH=${GENOME_DIR}/seqdata/dna.ext_1.0k.gff3
COUNTS_DIR=${BATCH_DIR}/counts



cd ${PROJECT_DIR}
mkdir -p ${BAM_DIR}
cd ${CLEAN_FASTQ_DIR}
for fastq_filepath in `ls *.clean.fastq.gz`; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`

    echo 'start mapping ...'
    echo `date`
    time ${BIN}/hisat2 -p ${nCPU} --no-spliced-alignment -x ${GENOME_INDEX} \
                       -U ${fastq_filepath} -S ${BAM_DIR}/${fastq_filename}.sam
    echo `date`
    echo 'finished.'

    cd ${BAM_DIR}
    ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
    ${BIN}/samtools index -c ${fastq_filename}.bam

    rm ${fastq_filename}.sam
    cd -
done


mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                     -o ${COUNTS_DIR}/cs.counts.gene.ext_1.0k.tsv 20181109*.bam




