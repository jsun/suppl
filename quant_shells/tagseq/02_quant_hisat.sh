#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_hisat
#$ -mods l_hard h_rt 288:00:00



pHISAT=1
pQuant=1
nCPU=16



PROJECT_DIR=~/projects/HuoberBrezel
BIN=~/local/bin
DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bam

GENOME_DIR=~/projects/db/genome/IWGSC_RefSeq_v2.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
GTF_FILEPATH=${GENOME_DIR}/seqdata/dna.gff3
COUNTS_DIR=${DATA_DIR}/counts



cd ${PROJECT_DIR}


if [ ${pHISAT} -eq 1 ]; then
# 
# mapping reads with HISAT2
# 
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
fi





if [ ${pQuant} -eq 1 ]; then
#
# count mapped reads
#
mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                     -o ${COUNTS_DIR}/tcs.counts.gene.tsv TaeRS2728_*.bam
${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                     -o ${COUNTS_DIR}/cs.counts.gene.tsv 20181109*.bam
gff_versions=("ext_0.5k" "ext_1.0k" "ext_1.5k" "ext_2.0k" "ext_2.5k" "ext_3.0k" "ext_3.5k" "ext_4.0k")
for gff in ${gff_versions[@]}; do
    ${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH%.gff3}.${gff}.gff3 -s 1 \
                         -o ${COUNTS_DIR}/tcs.counts.gene.${gff}.tsv TaeRS2728_*.bam
    ${BIN}/featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH%.gff3}.${gff}.gff3 -s 1 \
                         -o ${COUNTS_DIR}/cs.counts.gene.${gff}.tsv 20181109*.bam
done
fi


