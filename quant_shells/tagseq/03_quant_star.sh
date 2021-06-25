#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_star
#$ -mods l_hard h_rt 288:00:00


pSTAR=1
pQuant=1
nCPU=16



PROJECT_DIR=~/projects/HuoberBrezel
BIN=~/local/bin
DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bamstar


GENOME_DIR=~/projects/db/genome/IWGSC_RefSeq_v2.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_star

GTF_FILEPATH=${GENOME_DIR}/seqdata/dna.ext_1.0k.gff3
COUNTS_DIR=${DATA_DIR}/countsstar




cd ${PROJECT_DIR}


if [ ${pSTAR} -eq 1 ]; then
# 
# mapping reads with STAR
# 
mkdir -p ${BAM_DIR}
cd ${BAM_DIR}
for fastq_filepath in `ls ${CLEAN_FASTQ_DIR}/*.clean.fastq.gz`; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`

    echo 'start mapping ...'
    echo `date`

    time ${BIN}/STAR --runThreadN ${nCPU} \
                     --genomeDir ${GENOME_INDEX} \
                     --readFilesCommand zcat \
                     --readFilesIn ${fastq_filepath} \
                     --genomeLoad NoSharedMemory \
                     --alignIntronMax 1 \
                     --outMultimapperOrder Random \
                     --outFileNamePrefix ${BAM_DIR}/${fastq_filename}_
    
    echo `date`
    echo 'finished.'
    
    ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}_Aligned.out.sam
    ${BIN}/samtools index -c ${fastq_filename}.bam
    rm ${fastq_filename}_Aligned.out.sam
done
fi





if [ ${pQuant} -eq 1 ]; then
#
# count mapped reads
#
mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                     -o ${COUNTS_DIR}/tcs.counts.gene.ext_1.0k.tsv TaeRS2728_*.bam
${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                     -o ${COUNTS_DIR}/cs.counts.gene.ext_1.0k.tsv 20181109*.bam
fi


